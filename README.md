# Ranked-coal 
can be used to compute probabilities of ranked or unranked phylogenetic gene trees given a species tree under coalescent process.  

## Installation
After downloading the source code, type
```
make
```
This will create an executable called ranked, which can be futher run with some input options listed below.

## Usage
Program options:

| Option        | Input files   | Output files                   |
| ------------- |:-------------| :------------------------------|
| -rprob        | <ul><li>input a species tree file</li><li>input a file containing ranked gene trees (with branch lengths)</li><li> input a file containing gene tree topologies (optional)</li></ul>|<ul><li>STtopo.txt</li><li>probForEachGT.txt</li></ul>|
| -uprob        | <ul><li>input a species tree file</li><li>input a file containing unranked gene trees (without branch lengths)</li></ul>| <ul><li>STtopo.txt</li><li>outUnrGT.txt</li></ul>|

All input files should be in the Newick format. All trees are treated as rooted binary trees. We assume an ultrametric species tree (leaves of the tree are all equidistant from the root).  

When using *-rprob* option, the input may look like this:
```
./ranked -rprob <species-tree-file-name> <ranked-gene-tree-file-name> <gene-tree-topology-file-name>
```
where 
> ```<ranked-gene-tree-file-name>``` contains one or more gene trees with specified branch lengths in the Newick format. The taxon names of gene trees should match the taxon names of the corresponding species tree.   
> ```<gene-tree-topology-file-name>``` (optional input file) contains corresponding gene tree topologies (see Examples). 

When using *-uprob* option, the input may look like this:
```
./ranked -uprob <species-tree-file-name> <unranked-gene-tree-file-name>
```
> where ```<unranked-gene-tree-file-name>``` contains one unranked gene trees without branch lengths in the Newick format. The taxon names of gene trees should match the taxon names of the corresponding species tree. The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology.   






## Examples

