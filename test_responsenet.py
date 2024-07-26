import pytest
import argparse
import responsenet as rn

args = argparse.Namespace(edges_file='data/inputs/test-edges.txt', sources_file='data/inputs/test-sources.txt', targets_file='data/inputs/test-targets.txt', output='test-out', gamma=10, include_st=False, verbose=True, output_log=True)

## functions from main(). There has to be a better way to do this.
sources = rn.parse_nodes(args.sources_file)
targets = rn.parse_nodes(args.targets_file)

global _verbose 
global _include_st 
global _output_log
_verbose = args.verbose
_include_st = args.include_st
_output_log = args.output_log

gamma = args.gamma

G = rn.construct_digraph(args.edges_file)
G = rn.add_sources_and_targets(G, sources, targets)

out_file = args.output+"_gamma"+str(gamma)+".tsv"
out_log = args.output + out_file[6:-4] + ".log"

solver = rn.responsenet(G, gamma, out_file, out_log)

## also add test_lp() that counts the number of constraints
## optionally, add test functions that confirm that the sources/targest are being added appropriately.

def test_objective():
    assert round(solver.Objective().Value(),4) == -8.5729

def test_flows():
    for u,v in G.edges:
        if (u=="source" and v =="A") or (u=="A" and v=="B") or (u=="B" and v=="E") or (u=="E" and v=="target"):
            assert G[u][v]["flow"].solution_value() == 1.0
        else:
            G[u][v]["flow"].solution_value() == 0
