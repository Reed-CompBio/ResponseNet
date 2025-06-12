# ResponseNet approaches pathway reconstruction by modelling the problem as a minimum-cost flow optimization problem.
# However, unlike other ways of modelling this (e.g. PCSTs), ResponseNet optimizes flow by encoding the entire graph
# as an ILP problem.
#
# The algorithm for ResponseNet implemented below is described in the _Linear programming formulation_ section of the
# paper linked in `README.md`.
# 
# The genetic hits are the sources, and the differentially expressed genes are the targets. These node sets are
# associated with a weighed interactome, or a weighed, directed graph.

import argparse
import logging
import math
import networkx as nx
from ortools.linear_solver import pywraplp
from pathlib import Path
import warnings

# Global Variables that args can modify
_include_st = False
_output_log = False

def as_cost(weight: float) -> float:
    # We get the negated log of the weight as our cost: weight is from (0, 1] (truncated to 0.7)
    # where greater values are edges we want to keep in the interactome. This reframes
    # the 'weight maximization' game into a 'cost minimization' game. -log(...) will transform our weights
    # where lower is higher and higher is lower.
    return math.log(weight) * -1

def parse_nodes(node_file: Path):
    """ 
    Parse a list of sources or targets and return a set 
    
    @param node_file: the PATH file for a list of nodes
    @return: set of all nodes listed in file
    """
    lines = node_file.read_text().splitlines()
    nodes = set(map(str.strip, lines))
    return nodes

def construct_digraph(edges_file: Path, default_capacity=1):
    """
    Similar to MinCostFlow, we need to parse a list of undirected edges and 
    returns a graph object
    
    @param edges_file: the PATH file for an interactome
    @param default_capacity: the capacity (c_(ij)) to give to all of the edges initially.
    
    @return: the constructed graph object
    """
    
    ## Make a directed graph object.
    G = nx.DiGraph()
    
    # Go through edge_file, assign each node an id
    with open(edges_file) as edges_f:
        for line in edges_f:
            tokens = line.strip().split('\t')
            if len(tokens) != 3:
                raise ValueError(f"Provided line {line} does not have 3 tab-separated entries.")
            source, target, weight = tokens
            
            if not source in G:
                G.add_node(source)
            if not target in G:
                G.add_node(target)
           
            weight = float(weight)
            # As described in the paper (at "Weighting scheme for interactome edges"),
            # we truncate scores to be between 0 and 0.7: 
            #     "Because high edge weights could indicate unusually well-studied proteins or imperfectness 
            #     of the assumption of conditional independence, all weights were capped to a maximum value of 0.7"
            if weight > 0.7:
                weight = 0.7
            
            # Zero-weight or negative edges cause problems - note that we will take the negated log of the weight later.
            # TODO: can we do anything about zero-weight edges?
            if weight <= 0.0:
                warnings.warn(f"Edge {source} --> {target} has weight <= 0 ({weight}), this will cause problems.")
            
            cost = as_cost(weight)

            ## AR change "cost" to "weight" so it accurately reflects the value. 
            G.add_edge(source,
                        target,
                        cost=cost,
                        cap=default_capacity)

    return G
    
def add_sources_and_targets(G: nx.DiGraph, sources: set[str], targets: set[str]) -> nx.DiGraph:
    """
    Add a 'false' super source and target node to the DiGraph, helpful
    for organization and essential to the ILP.

    @param G: DiGraph object
    @param sources: set of all source nodes
    @param targets: set of all target nodes
        
    @return: modified DiGraph object with faux source and target
    """

    # Divide the capacity evently across the sources and targets.
    source_weight = 1 / len(sources)
    target_weight = 1 / len(targets)
    
    source_cap = source_weight
    target_cap = target_weight

    source_cost = as_cost(source_weight)
    target_cost = as_cost(target_weight)

    if G.has_node("source"):
        raise ValueError("A node named 'source' is already present - ResponseNet can't add a super-source node.")
    if G.has_node("target"):
        raise ValueError("A node named 'target' is already present - ResponseNet can't add a super-target node.")
    
    G.add_node("source")
    G.add_node("target")

    for source in sources:
        logging.debug(f'Looping through source: {source}')
        if source in G:
            G.add_edge("source",
                        source,
                        cost=source_cost,
                        cap=source_cap)
        else:
            warnings.warn(f"Source '{source}' not found in graph")

    for target in targets:
        logging.debug(f'Looping through target: {target}')
        if target in G:
            G.add_edge(target,
                        "target",
                        cost=target_cost,
                        cap=target_cap)
        else:
            warnings.warn(f"Target '{target}' not found in graph")   
            
    return G
    
def prepare_variables(solver: pywraplp.Solver, G: nx.DiGraph) -> dict[tuple, pywraplp.Variable]:
    """
    This section systematically creates variables for the ILP and saves them
    both in a dictionary and as an attribute for each edge in G

    @param solver: solver object that the LP depends on
    @param G: graph object of interactome
        
    @returns flows: dictionary of all variables in the solver
    """

    # Here, we want to construct all of the flow variables for the ILP.
    # We take all of the present "cap" values per edge (see `construct_digraph`)
    # and add them to the ILP to be optimized.
    # This adds the c_(ij) step present in the paper,
    # but makes them flexible variables who can be between 0 and cap.
    flows: dict[tuple, pywraplp.Variable] = dict()
    extras = 0
    for i, j in G.edges():
        edge = (i,j)
        if edge not in flows:
            flows[edge] = solver.NumVar(0.0, G[i][j]["cap"], f"Flows{edge}")
            G.get_edge_data(i, j)["flow"] = flows[edge]
        else:
            logging.info("Found repeating edge: {edge}")
            extras += 1
    logging.info(f"There were {extras} repeat edges.")

    # [On debug mode] log the status of the solver
    debug_log_solver(solver)
    
    return flows
    
def prepare_constraints(solver: pywraplp.Solver, G: nx.DiGraph) -> list[pywraplp.Constraint]:
    """
    This section systematically applies constraints on each node and all edges
    to make sure that any flow entering a node also exits a node

    @param solver: solver object that LP depends on
    @param G: graph object of interactome

    @return constraints: list object containing all constraints in the LP
    """

    constraints: list[pywraplp.Constraint] = []
    for i, node in enumerate(G.nodes):
        if node == "source" or node == "target":
            continue

        # Creating constraint for each node, named after the node.
        # We establish that the node must be constrained to the value zero, as to
        # say that all of the coefficients attached to the constraint must add up to zero.
        curr_constraint = solver.Constraint(0.0, 0.0, node)
        
        constraints.append(curr_constraint)
        G.nodes[node]["constraint"] = curr_constraint

        in_edges = G.in_edges(node)
        out_edges = G.out_edges(node)

        # Since the node must have a final value of zero,
        # we add 1 and -1 coefficients to the incoming and outgoing edges, respectively,
        # to say that, as the top-level docstring implies, all flow entering this node
        # also exits it.
        for u,v in in_edges:
            assert v == node
            constraints[i].SetCoefficient(G[u][v]["flow"], 1)
            
        for u,v in out_edges:
            assert u == node
            constraints[i].SetCoefficient(G[u][v]["flow"], -1)

    # Adding a final constraint to make sure all flows going from the source
    # and to the target are equivalent. The same idea for this constraint
    # is present in the above for-loop.
    constraints.append(solver.Constraint(0.0, 0.0, "source"))

    for j, k in G.out_edges("source"):
        constraints[-1].SetCoefficient(G[j][k]["flow"], 1)
    for j, k in G.in_edges("target"):
        constraints[-1].SetCoefficient(G[j][k]["flow"], -1)
        
    # [On debug mode] log the status of the solver
    debug_log_solver(solver)

    return constraints
            
def prepare_objective(solver: pywraplp.Solver, G: nx.DiGraph, flows: dict, gamma: int) -> pywraplp.Objective:
    """
    This segment goes through all edges in the graph and sets a coefficient on each variable in the LP

    @param solver: solver object that LP depends on
    @param G: graph object of interactome
    @param flows: dictionary of all flow variables for the solver
    @param gamma: user defined value that determines graph size
    
    @returns objective: solver objective with all constraints
    """
    objective: pywraplp.Objective = solver.Objective()
    
    # The general goal of this objective is to minimize
    # the cost of all of the edges.
    for i,j in G.edges():
        # We want to minimize our costs (see `as_cost` for how weights are processed into costs)
        cost = G[i][j]["cost"] 
        if i == "source":
            # The higher gamma is, the more flow that is allowed to transfer through the network,
            # as this subtraction rewards any flow that travels from the super-source through the
            # rest of the network.
            cost = cost - gamma
        objective.SetCoefficient(flows[i,j], cost)
    
    objective.SetMinimization()
    
    # [On debug mode] log the status of the solver
    debug_log_solver(solver)

    return objective  

def debug_log_solver(solver):
    """
    Helper function to print contents of solver (constraints, variables, objective) for debugging
    """

    logging.debug('**' * 25)
    logging.debug(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','))
    logging.debug('**' * 25)

## AR make this return the solver, for testing.
def responsenet(G: nx.DiGraph, gamma: int, out_file: Path, out_log: Path) -> pywraplp.Solver:
    """ 
    The NEW ILP solver for ResponseNet, using GLOP.

    @param G: graph object of interactome
    @param gamma: user defined integer determining size of output graph
    @param out_file: path to the output file for writing the LP solution
    """
    
    solver: pywraplp.Solver = pywraplp.Solver.CreateSolver("GLOP")
    if not solver:
        raise RuntimeError("Could not construct GLOP solver.")
        
    # Data structures that define the ILP, kept for your debugging pleasure
    flows = prepare_variables(solver, G)
    _constraints = prepare_constraints(solver, G)
    _objective = prepare_objective(solver, G, flows, gamma)
    
    print("Attempting solve of flows...")
    status = solver.Solve()
    
    if status == pywraplp.Solver.OPTIMAL:
        print("Solved! \n")
    else:
        raise RuntimeError("The problem does not have an optimal solution.")
    
    write_output_to_tsv(G, solver, out_file, out_log)
    return solver

def write_output_to_tsv(G: nx.DiGraph, solver: pywraplp.Solver, out_file: Path, out_log: Path):
    '''
    Write output of solver.Solve() over graph obj to an output file specified 
    by out_file
    
    @G : graph object
    @solver: contains the answer to the LP
    @out_file: Path to output file
    @out_log: Path to output log
    '''
    with out_file.open("w") as output_f:
        print(f"Objective value = {solver.Objective().Value():0.1f}")
        print(f"Solved in {(float(solver.wall_time())/1000)} seconds")
        
        output_f.write("Interactor 1" + '\t' + "Interactor 2" + '\t' + "Flow" + "\n")
        for u,v in G.edges:
            # Check to see if we want to actually include the artificial source and target  
            if (u == "source" or v == "target") and not _include_st:
                continue
            else:
                if G[u][v]["flow"].solution_value() > 0.0 and G[u][v]["flow"].solution_value() <= 1.0:   
                    output_f.write(str(u)+"\t"+str(v)+"\t"+str(G[u][v]["flow"].solution_value())+"\n")

    # Format for output log, including the entire solver information
    if _output_log:
        with out_log.open("w") as out_l:
            out_l.write("Objective value = " + str(solver.Objective().Value()) + '\n')
            out_l.write("Solved in " + str(float(solver.wall_time()) / 1000) + " seconds" + '\n\n')
            out_l.write("Solver:\n")
            out_l.write(str(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ',')))

def main(args):
    print("Running ResponseNet...")

    sources = parse_nodes(Path(args.sources_file))
    targets = parse_nodes(Path(args.targets_file))
    
    # Modifying global variables based on args
    global _include_st 
    global _output_log
    _include_st = args.include_st
    _output_log = args.output_log

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    gamma = args.gamma
    
    G = construct_digraph(args.edges_file)
    
    G = add_sources_and_targets(G, sources, targets)
    
    # AR make this a TXT file. Keep the same formatting. Should we have headers
    out_file = args.output + "_gamma" + str(gamma) + ".txt"
    out_log = args.output  + "_gamma" + str(gamma) + ".log"
    responsenet(G, gamma, Path(out_file), Path(out_log))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--edges_file',
                        help='Network file. File should be in SIF file format.',
                        type=str,
                        required=True)
    parser.add_argument('--sources_file',
                        help='File which denotes source nodes, with one node per line.',
                        type=str,
                        required=True)
    parser.add_argument('--targets_file',
                        help='File which denotes source nodes, with one node per line.',
                        type=str,
                        required=True)
    parser.add_argument('--output',
                        help='Prefix for all output files.',
                        type=str,
                        required=True)
    parser.add_argument('--gamma',
                        help='The size of the output graph. Default = 10.',
                        type=int,
                        required=False,
                        default=10)
    parser.add_argument('-st','--include_st',
                        help='Determines whether output should include artificial Source and Target nodes. By default does not include them.',
                        action='store_true')
    parser.add_argument('-v','--verbose',
                        help='Include verbose console output',
                        action='store_true')
    parser.add_argument('-o', '--output_log',
                        help='Create output log',
                        action='store_true')

    args = parser.parse_args()
    print(args)

    main(args)
