# ResponseNet
Implementation of the ResponseNet algorithm. The algorithm can be found in this paper:

[Bridging high-throughput genetic and transcriptional data reveals cellular responses to alpha-synuclein toxicity](https://www.nature.com/articles/ng.337) by Yeger-Lotem et al., Nature Genetics 2009.

This method has appeared in various webservers (most recently as [ResponseNet v3](https://pubmed.ncbi.nlm.nih.gov/31114913/)), but the underlying algorithm is not publicly available. 

### Usage 
```
usage: responsenet.py [-h] --edges_file EDGES_FILE --sources_file SOURCES_FILE --targets_file TARGETS_FILE --output OUTPUT
                      [--gamma GAMMA] [-st] [-v] [-o]

options:
  -h, --help            show this help message and exit
  --edges_file EDGES_FILE
                        Network file. File should be in SIF file format.
  --sources_file SOURCES_FILE
                        File which denotes source nodes, with one node per line.
  --targets_file TARGETS_FILE
                        File which denotes source nodes, with one node per line.
  --output OUTPUT       Prefix for all output files.
  --gamma GAMMA         The size of the output graph. Default = 10
  -st, --include_st     Determines whether output should include artificial Source and Target nodes
  -v, --verbose         Include verbose console output
  -o, --output_log      Create output log
  ```

### Quick Start

We have provided a small example of a directed weighted graph, with a source `A` and a target `E`:

![JPEG image-4B8D-A0E6-92-0](https://github.com/user-attachments/assets/51260387-bf38-4659-9577-ff1484937af0)

You can run the example with the command:

```
python responsenet.py --edges_file data/inputs/test-edges.txt --sources_file data/inputs/test-sources.txt --targets_file data/inputs/test-targets.txt --output out -v -o
```

This will output the final file of edges with non-zero flow as well as a log file (`-o`) with the LP. It will also output some information to the console (`-v`). The ouput edges file will look like

```
A	B	1.0
B	E	1.0
```

Which corresponds to these edges in the graph:

![IMG_F2C39098471A-1](https://github.com/user-attachments/assets/16fa0904-c048-4b68-a429-f4b7df498f26)

### Some notes:
- Edge weights are positive are should range between 0 and 1. Similar to the paper, edge weights greater than 0.7 are truncated to 0.7.
- The gamma variable determines the size of the network. The default setting is 10, same as the paper.
