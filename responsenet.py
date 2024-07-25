# The beginning of the rewrite of responsenet

from ortools.linear_solver import pywraplp
import argparse
import networkx as nx
import math

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
    returns a graph object and idDict
    
    Parameters:
        @edges_file : PATH()
            the PATH file for an interactome
    
    Returns:
        @G: graph object
        @idDict: Dictionary of all nodes in interactome, 
        mapped to an integer value for MCF
    """
    
    #G = min_cost_flow.SimpleMinCostFlow()
    G = nx.DiGraph()
    idDict = dict()
    curID = 0

    weights = []
    all_ones = 0
    
    # Go through edge_file, assign each node an id
    with open(edges_file) as edges_f:
        for line in edges_f:
            tokens = line.strip().split()
            node1 = tokens[0]
            
            # Add nodes to idDict if they aren't there
            # NOTE: This is an artifact of MCF, idDict functionality may be modified or removed
            if not node1 in idDict:
                idDict[node1] = curID
                curID += 1
                G.add_node(node1, ident = curID)
            node2 = tokens[1]
            if not node2 in idDict:
                idDict[node2] = curID
                curID += 1
                G.add_node(node2, ident = curID)
           
        # TODO: Add check for all weights to be between 0,1
        #       Do an error if all weights are 1
        #       NOTE: Check *should* be done, will test later
            w = float(tokens[2])
            if w <= 1.0 and w >= 0.0:
                weights.append(w)
                if w == 1.0:
                    all_ones += 1
            else:
                raise(ValueError)
            
            # NOTE: Commenting out as part of idDict depreciation, will remove later
            # G.add_edge(idDict[node1], 
            #            idDict[node2], 
            #            cost = w, 
            #            cap = default_capacity)

            G.add_edge(node1,
                        node2,
                        cost = w,
                        cap = default_capacity)
        
        # Will throw error if all weights are equal to 1 or 1.0
        if all_ones == len(weights):
            raise(ValueError)

        idDict["maxID"] = curID
        return G, idDict
    
def add_sources_and_targets(G, sources, targets, idDict):
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
        @idDict : dict()
            dictionary that assigns integer values to each node
            to uniquely identify each node

    Returns:
        @G : modified DiGraph object with faux source and target
    """
    source_weight = 1/len(sources)
    target_weight = 1/len(targets)
    
    source_cap = source_weight
    target_cap = target_weight
    
    curID = idDict["maxID"]
    idDict["source"] = curID
    G.add_node("source", ident = curID)
    curID += 1
    idDict["target"] = curID
    G.add_node("target", ident = curID)

    for source in sources:
        print(source)
        if source in G:
            print("source found")
            # Commenting out to also prepare for loss of idDict
            # G.add_edge(idDict["source"], 
            #            idDict[source], 
            #            cost = source_weight, 
            #            cap = source_cap)
            
            G.add_edge("source",
                        source,
                        cost = source_weight,
                        cap = source_cap)

        else:
            print(f"Source: {source} not found in graph")

    for target in targets:
        print(target)
        if target in idDict:
            print("target found")
            # See above comment
            # G.add_edge(idDict[target], 
            #            idDict["target"], 
            #            cost = target_weight, 
            #            cap = target_cap)
           
            G.add_edge(target,
                        "target",
                        cost = target_weight,
                        cap = target_cap)

        else:
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
            # Seeing how this changes things (it does)
            flows[edge] = solver.NumVar(0.0, G[i][j]["cap"], f"Flows{edge}")
            G.get_edge_data(edge[0],edge[1])["flow"] = flows[edge]
        else:
            print("repeat")
            print(edge)
            extras += 1
    print(f"We had {extras} repeat edges")

    # Helpful debugging statement for LP solver    
    # print('**'*25)
    # print(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','), sep='\n')
    # print('**'*25)
    
    return flows
    
def prepare_constraints(solver, G, idDict):
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
        
        if G.nodes[node]["ident"] == idDict["source"] or G.nodes[node]["ident"] == idDict["target"]:
            continue   
        
        # 
        assert(i == idDict[node])

        # Trying out a new way of marking constraints, while also wrapping the data into the G object
        #curr_constraint = solver.Constraint(idDict[node],solver.infinity())
        # Trying to figure out if line 228 and line 230 make any differences
        curr_constraint = solver.Constraint(i, solver.infinity())
        
        constraints.append(curr_constraint)
        G.nodes[node]["constraint"] = curr_constraint

        for u,v in in_edges:
            assert v == node
            constraints[i].SetCoefficient(G[u][v]["flow"],1)
            
        for u,v in out_edges:
            assert u == node
            constraints[i].SetCoefficient(G[u][v]["flow"],-1)
            
        constraints[i].SetBounds(0,0)


    # Modified for depreciation of idDict
    constraints.append(solver.Constraint(idDict["source"], solver.infinity()))

    # Testing if we should 
    
    # trying something silly
    # constraints.append(solver.Constraint(0.0, solver.infinity()))

    for j,k in list(G.out_edges("source")):
        constraints[-1].SetCoefficient(G[j][k]["flow"],1)
    for j,k in list(G.in_edges("target")):
        constraints[-1].SetCoefficient(G[j][k]["flow"],-1)
        
    constraints[-1].SetBounds(0,0)

    # Helpful debugging statement for LP    
    # print('**'*25)
    # print(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','), sep='\n')
    # print('**'*25)
    
    return constraints
            
def prepare_objective(solver, G, flows, gamma, s):
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
        @s : int()
            integer stand-in for the "source" node
    
    Returns:
        @objective : solver objective with all constraints
    """
    objective = solver.Objective()
    
    for i,j in G.edges():
        
        log_weight = (math.log(G[i][j]["cost"])) * (-1)
        
        if i == "source":
            log_weight = log_weight - gamma  
            print("adjusting for source")
        objective.SetCoefficient(flows[i,j], log_weight) 
    
    objective.SetMinimization()
        
    # Helpful debugging statement to show status of LP solver
    print('**'*25)
    print(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','), sep='\n')
    print('**'*25)

    return objective  


def responsenet(G, idDict, gamma, out_file):
    """ 
    The NEW ILP solver for ResponseNet, using GLOP.

    Parameters:
        @G : nx.DiGraph()
            graph object of interactome
        @idDict : dict()
            dictionary of all nodes in network
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
    
    s = idDict["source"]
    
    # Data structures that define the ILP, kept for your debugging pleasure
    flows = prepare_variables(solver, G)
    constraints = prepare_constraints(solver, G, idDict)
    objective = prepare_objective(solver, G, flows, gamma, s)
    
    print("Attempting solve of flows")
    status = solver.Solve()
    
    if status == pywraplp.Solver.OPTIMAL:
        print("Solved!")

    else:
        print("The problem does not have an optimal solution.")
        return
    
    write_output_to_tsv(status, G, solver, out_file)
    
def write_output_to_tsv(status, G, solver, out_file, include_st = False):
    '''
    Write output of solver.Solve() over graph obj to an output file specified 
    by out_file
    
    Params:
        @status: status of solver.Solve()

    # TODO: add check for if include_st == True/False
    # TODO: Confirm if SIF is preferred output type
    '''
    with open(out_file, "w") as output_f:
        print(f"Objective value = {solver.Objective().Value():0.1f}")
        
        for u,v in G.edges:    
            if G[u][v]["flow"].solution_value() > 0.0 and G[u][v]["flow"].solution_value() < 1.0:   
                #print(G[u][v]["flow"],'-->',G[u][v]["flow"].solution_value())
                output_f.write(str(u)+"\t"+str(v)+"\t"+str(G[u][v]["flow"].solution_value())+"\n")
    return

def main(args):
    
    sources = parse_nodes(args.sources_file)
    targets = parse_nodes(args.targets_file)
    
    gamma = args.gamma
    
    G, idDict = construct_digraph(args.edges_file)
    
    G = add_sources_and_targets(G, sources, targets, idDict)
    
    out_file = args.output+"_gamma"+str(gamma)+".tsv"

    responsenet(G, idDict, gamma, out_file)
    
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
                        help='The size of the output graph',
                        type=int,
                        required=True,
                        default=10)

    args = parser.parse_args()
main(args)
