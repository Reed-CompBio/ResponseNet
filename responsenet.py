# The beginning of the rewrite of responsenet

from ortools.linear_solver import pywraplp
import argparse
import networkx as nx
import math
import warnings

# Global Variables that args can modify
_verbose = False
_include_st = False
_output_log = False

def parse_nodes(node_file):
    """ 
    Parse a list of sources or targets and return a set 
    
    Parameters:
        @node_file : PATH()
            the PATH file for a list of nodes

    Returns:
        @nodes: set of all nodes listed in file
    """
    with open(node_file) as node_f:
        lines = node_f.readlines()
        nodes = set(map(str.strip, lines))
    return nodes

def construct_digraph(edges_file, default_capacity= 1):
    """
    Similar to MinCostFlow, we need to parse a list of undirected edges and 
    returns a graph object
    
    Parameters:
        @edges_file : PATH()
            the PATH file for an interactome
    
    Returns:
        @G: graph object
    """
    
    ## Make a directed graph object.
    G = nx.DiGraph()
    
    # Go through edge_file, assign each node an id
    with open(edges_file) as edges_f:
        for line in edges_f:
            tokens = line.strip().split()
            node1 = tokens[0]
            
            if not node1 in G:
                G.add_node(node1)
            node2 = tokens[1]
            if not node2 in G:
                G.add_node(node2)
           
            w = float(tokens[2])
            # From the paper: truncate scores to be between 0 and 0.7. 
            # Because high edge weights could indicate unusually well-studied proteins or imperfectness 
            # of the assumption of conditional independence, all weights were capped to a maximum value of 0.7"
            if w > 0.7:
                w = 0.7
            
            # zero-weight or negative edges will cause a problem. 
            if w <= 0.0:
                warnings.warn(f"Edge {tokens[0]} --> {tokens[1]} has weight <= 0, this will cause problems")

            ## AR change "cost" to "weight" so it accurately reflects the value. 
            G.add_edge(node1,
                        node2,
                        cost = w,
                        cap = default_capacity)

        return G
    
def add_sources_and_targets(G, sources, targets):
    """
    Add a false source and target node to the DiGraph, helpful
    for organization and essential to the ILP

    Parameters:
        @G : nx.DiGraph()
            DiGraph object
        @sources : set()
            set of all source nodes
        @targets : set()
            set of all target nodes
        
    Returns:
        @G : modified DiGraph object with faux source and target
    """

    # Divide the capacity evently across the sources and targets
    source_weight = 1/len(sources)
    target_weight = 1/len(targets)
    
    source_cap = source_weight
    target_cap = target_weight
    
    G.add_node("source")
    G.add_node("target")

    for source in sources:
        if _verbose:
            print(source)
        if source in G:
            if _verbose:
                print("source found")
            G.add_edge("source",
                        source,
                        cost = source_weight,
                        cap = source_cap)

        else:
            if _verbose:
                print(f"Source: {source} not found in graph")

    for target in targets:
        
        if _verbose:
            print(target)
        if target in G:
            if _verbose:
                print("target found")
            G.add_edge(target,
                        "target",
                        cost = target_weight,
                        cap = target_cap)

        else:
            if _verbose:
                print(f"Target: {target} not found in graph")   
            
    return G
    
def prepare_variables(solver, G):
    """
    This section systematically creates variables for the ILP and saves them
    both in a dictionary and as an attribute for each edge in G

    Parameters:
        @solver : pywraplp.Solver()
            solver object that the LP depends on
        @G : nx.DiGraph()
            graph object of interactome
        
    Returns:
        @flows: dictionary of all variables in the solver
    """
    flows = dict()
    extras = 0
    for i,j in G.edges():
        edge = (i,j)
        if edge not in flows:
            # Need to set max value for each edge to be the max capacity of given edge
            flows[edge] = solver.NumVar(0.0, G[i][j]["cap"], f"Flows{edge}")
            G.get_edge_data(edge[0],edge[1])["flow"] = flows[edge]
        else:
            if _verbose:
                print("repeat")
                print(edge)
            extras += 1
    if _verbose:
        print(f"We had {extras} repeat edges")

    # Helpful debugging statement for LP solver    
    # print_solver(solver)
    
    return flows
    
def prepare_constraints(solver, G):
    """
    This section systematically applies constraints on each node and all edges
    to make sure that any flow entering a node also exits a node

    Parameters:
        @solver : pywraplp.Solver()
            solver object that LP depends on
        @G : nx.DiGraph()
            graph object of interactome
        @idDict : dict()
            dictionary of all nodes in network

    Returns:
        @constraints: list object containing all constraints in the LP
    """
    constraints = []
    for i,  node in enumerate(G.nodes):
        
        in_edges = list(G.in_edges(node))
        out_edges = list(G.out_edges(node))
        
        if node == "source" or node == "target":
            continue   
        
        # Creating constraint for each node, constraint has bounds 0,0 
        # and is named after the node
        curr_constraint = solver.Constraint(0.0, 0.0, node)
        
        constraints.append(curr_constraint)
        G.nodes[node]["constraint"] = curr_constraint

        for u,v in in_edges:
            assert v == node
            constraints[i].SetCoefficient(G[u][v]["flow"],1)
            
        for u,v in out_edges:
            assert u == node
            constraints[i].SetCoefficient(G[u][v]["flow"],-1)

    # Adding a final constraint to make sure all flows going from the source
    # and to the target are equivalent
    constraints.append(solver.Constraint(0.0, 0.0, "source"))

    for j,k in list(G.out_edges("source")):
        constraints[-1].SetCoefficient(G[j][k]["flow"],1)
    for j,k in list(G.in_edges("target")):
        constraints[-1].SetCoefficient(G[j][k]["flow"],-1)
        
    # Helpful debugging statement for LP    
    # print_solver(solver)

    return constraints
            
def prepare_objective(solver, G, flows, gamma):
    """
    This segment goes through all edges in the graph and sets a coefficient on each variable in the LP

    Parameters:
        @solver : pywraplp.Solver()
            solver object that LP depends on
        @G : nx.DiGraph()
            graph object of interactome
        @flows : dict()
            dictionary of all flow variables for the solver
        @gamma : int()
            user defined value that determines graph size
    
    Returns:
        @objective : solver objective with all constraints
    """
    objective = solver.Objective()
    
    for i,j in G.edges():
        
        log_weight = (math.log(G[i][j]["cost"])) * (-1)
        
        if i == "source":
            log_weight = log_weight - gamma  
            if _verbose:
                print("adjusting for source")
        objective.SetCoefficient(flows[i,j], log_weight) 
    
    objective.SetMinimization()
    
    # Helpful debugging statement to show status of LP solver
    # print_solver(solver)

    return objective  

def print_solver(solver):
    """
    Helper function to print contents of solver (constraints, variables, objective) for debugging
    """

    print('**'*25)
    print(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','), sep='\n')
    print('**'*25)

## AR make this return the solver, for testing.
def responsenet(G, gamma, out_file, out_log):
    """ 
    The NEW ILP solver for ResponseNet, using GLOP.

    Parameters:
        @G : nx.DiGraph()
            graph object of interactome
        @gamma : int()
            user defined integer determining size of output graph
        @out_file : PATH()
            PATH to the output file for writing the LP solution
    
    Returns:
        Nothing
    """
    
    solver = pywraplp.Solver.CreateSolver("GLOP")
    if not solver:
        return
        
    # Data structures that define the ILP, kept for your debugging pleasure
    flows = prepare_variables(solver, G)
    constraints = prepare_constraints(solver, G)
    objective = prepare_objective(solver, G, flows, gamma)
    
    print("Attempting solve of flows...")
    status = solver.Solve()
    
    if status == pywraplp.Solver.OPTIMAL:
        print("Solved! \n")

    else:
        print("The problem does not have an optimal solution.")
        return
    
    write_output_to_tsv(G, solver, out_file, out_log)
    return solver

def write_output_to_tsv(G, solver, out_file, out_log):
    '''
    Write output of solver.Solve() over graph obj to an output file specified 
    by out_file
    
    Params:
        @G : graph object
        @solver : pywraplp.Solver() object, contains the answer to the LP
        @out_file : str of output file name/path
        @out_log : str of output log file name/path
    '''
    with open(out_file, "w") as output_f:
        print(f"Objective value = {solver.Objective().Value():0.1f}")
        print(f"Solved in {(float(solver.wall_time())/1000)} seconds")
        
        output_f.write("Interactor 1" + '\t' + "Interactor 2" + '\t' + "Flow" + '\n')
        for u,v in G.edges:

            # Check to see if we want to actually include the artificial source and target  
            if (u == "source" or v == "target") and not _include_st:
                continue
            else:
                if G[u][v]["flow"].solution_value() > 0.0 and G[u][v]["flow"].solution_value() <= 1.0:   
                    output_f.write(str(u)+"\t"+str(v)+"\t"+str(G[u][v]["flow"].solution_value())+"\n")

    # Format for output log, including the entire solver information
    if _output_log:
        with open(out_log, "w") as out_l:
            out_l.write("Objective value = " + str(solver.Objective().Value()) +'\n')
            out_l.write("Solved in " + str(float(solver.wall_time())/1000) + " seconds"+'\n\n')
            out_l.write("Solver:\n")
            out_l.write(str(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ',')))

        return

def main(args):
    
    print("Running ResponseNet...")

    sources = parse_nodes(args.sources_file)
    targets = parse_nodes(args.targets_file)
    
    # Modifying global variables based on args
    global _verbose 
    global _include_st 
    global _output_log
    _verbose = args.verbose
    _include_st = args.include_st
    _output_log = args.output_log

    gamma = args.gamma
    
    G = construct_digraph(args.edges_file)
    
    G = add_sources_and_targets(G, sources, targets)
    
    # AR make this a TXT file. Keep the same formatting. Should we have headers
    out_file = args.output+"_gamma"+str(gamma)+".txt"
    out_log = args.output +"_gamma"+str(gamma) + ".log"
    responsenet(G, gamma, out_file, out_log)
    
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
