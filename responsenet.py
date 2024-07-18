# The beginning of the rewrite of responsenet

from ortools.linear_solver import pywraplp
import argparse
import networkx as nx
import math

def parse_nodes(node_file):
    ''' Parse a list of sources or targets and return a set '''
    with open(node_file) as node_f:
        lines = node_f.readlines()
        nodes = set(map(str.strip, lines))
    return nodes

def construct_digraph(edges_file):
    """
    Similar to MinCostFlow, we need to parse a list of undirected edges and returns a graph object and idDict
    
    Parameters:
        edges_file : PATH()
            the PATH file for an interactome
    
    Returns:
        @G: graph object
        @idDict: Dictionary of all nodes in interactome, mapped to an integer value for MCF
    """
    
    #G = min_cost_flow.SimpleMinCostFlow()
    G = nx.DiGraph()
    idDict = dict()
    curID = 0
    default_capacity = 1
    
    # Go through edge_file, assign each node an id
    with open(edges_file) as edges_f:
        for line in edges_f:
            tokens = line.strip().split()
            node1 = tokens[0]
            
            # Add nodes to idDict if they aren't there
            if not node1 in idDict:
                idDict[node1] = curID
                curID += 1
            node2 = tokens[1]
            if not node2 in idDict:
                idDict[node2] = curID
                curID += 1
           
        # TODO: Add check for all weights to be between 0,1
        #       Do an error if all weights are 1
            w = float(tokens[2])
            
            G.add_edge(idDict[node1], idDict[node2], cost = w, cap = default_capacity)
            
        idDict["maxID"] = curID
        return G, idDict
    
def add_sources_and_targets(G, sources, targets, idDict, flow):
    """
    """
    source_weight = 1/len(sources)
    target_weight = 1/len(targets)
    
    G1 = G.copy()
    
    source_cap = source_weight
    target_cap = target_weight
    
    # subsets capturing the source and target nodes
    gen = []
    tra = []
    
    curID = idDict["maxID"]
    idDict["source"] = curID
    curID += 1
    idDict["target"] = curID

    for source in sources:
        print(source)
        if source in idDict:
            print("found")
            G1.add_edge(idDict["source"], idDict[source], cost = source_weight, cap = source_cap)
            gen.append(idDict[source])

    for target in targets:
        print(target)
        if target in idDict:
            G1.add_edge(idDict[target], idDict["target"], cost = target_weight, cap = target_cap)
            tra.append(idDict[target])
            
    return G1
    
def prepare_variables(solver, G1, default_capacity):
    flows = dict()
    extras = 0
    for edge in G1.edges():
        if edge not in flows:
            flows[edge] = solver.NumVar(0.0, default_capacity, f"Flows{edge}")
            G1.get_edge_data(edge[0],edge[1])["flow"] = flows[edge]
        else:
            print("repeat")
            print(edge)
            extras += 1
    print(f"We had {extras} repeat edges")
    
    print('**'*25)
    print(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','), sep='\n')
    print('**'*25)
    
    return flows
    
def prepare_constraints(solver, G1, idDict):
    constraints = []
    for i,  node in enumerate(G1.nodes):
        
        in_edges = list(G1.in_edges(node))
        out_edges = list(G1.out_edges(node))
        
        if node == idDict["source"] or node == idDict["target"]:
            continue
        
        constraints.append(solver.Constraint(node,solver.infinity()))
        
        for u,v in in_edges:
            assert v == node
            constraints[i].SetCoefficient(G1[u][v]["flow"],1)
            
        for u,v in out_edges:
            assert u == node
            constraints[i].SetCoefficient(G1[u][v]["flow"],-1)
            
        constraints[i].SetBounds(0,0)
    
    constraints.append(solver.Constraint(idDict["source"], solver.infinity()))
    
    for j,k in list(G1.out_edges(idDict["source"])):
        constraints[-1].SetCoefficient(G1[j][k]["flow"],1)
    for j,k in list(G1.in_edges(idDict["target"])):
        constraints[-1].SetCoefficient(G1[j][k]["flow"],-1)
        
    constraints[-1].SetBounds(0,0)

    
    print('**'*25)
    print(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','), sep='\n')
    print('**'*25)
    
    return constraints
            
def prepare_objective(solver, G1, flows, gamma, s):
    objective = solver.Objective()
    
    for i,j in G1.edges():
        
        log_weight = (math.log(G1[i][j]["cost"])) * (-1)
        
        if i == s:
            log_weight = log_weight - gamma  
            print("adjusting for source")
        objective.SetCoefficient(flows[i,j], log_weight) 
    
    objective.SetMinimization()
    
    return objective  
    
    print('**'*25)
    print(solver.ExportModelAsLpFormat(False).replace('\\', '').replace(',_', ','), sep='\n')
    print('**'*25)

def responsenet(G, G1, idDict, gamma, return_solver= False):
    """ The NEW ILP solver for MinCostFlow, using glop.
    """
    
    # temporary value, will make it changable later
    default_capacity = 1

    solver = pywraplp.Solver.CreateSolver("GLOP")
    if not solver:
        return
    
    s = idDict["source"]
    
    flows = prepare_variables(solver, G1, default_capacity)
    constraints = prepare_constraints(solver, G1, idDict)
    objective = prepare_objective(solver, G1, flows, gamma, s)
    
    if return_solver:
            return flows, constraints, objective, solver
    else:
        return None, None, None, solver
    
def write_output(status, G1, solver):
    if status == pywraplp.Solver.OPTIMAL:
        print("Solution:")
        print(f"Objective value = {solver.Objective().Value():0.1f}")
        for u,v in G1.edges:
            if G1[u][v]["flow"].solution_value() > 0.0:
                print(G1[u][v]["flow"],'-->',G1[u][v]["flow"].solution_value())
    else:
        print("The problem does not have an optimal solution.")    

def main():
    
    sources = parse_nodes("sources.txt")
    targets = parse_nodes("targets.txt")
    
    gamma = 10
    
    G, idDict = construct_digraph("edges.txt")
    
    G1 = add_sources_and_targets(G, sources, targets, idDict, 1)
    
    flows, constraints, objective, solver = responsenet(G, G1, idDict, gamma, True)
    
    return G, G1, idDict, sources, targets, solver, flows, constraints, objective
    
    
    