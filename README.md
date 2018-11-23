# Ranked-coal 
can be used to compute probabilities of ranked or unranked phylogenetic gene tree topologies given a species tree under coalescent process.  

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
| -rprob        | <ul><li>input a species tree file</li><li>input a file containing ranked gene trees (with branch lengths)</li><li> input a file containing gene tree topologies (optional)</li></ul>|<ul><li>STtopo.txt</li><li>outRankGT.txt</li></ul>|
| -uprob        | <ul><li>input a species tree file</li><li>input a file containing unranked gene trees (without branch lengths)</li></ul>| <ul><li>STtopo.txt</li><li>outUnrGT.txt</li></ul>|

All input files should be in the Newick format. All trees are treated as rooted binary trees. We assume an ultrametric species tree (leaves of the tree are all equidistant from the root). User can run *ranked-coal* with either *-rprob* or *-uprob* option as shown below.  

```
./ranked -rprob <species-tree-file-name> <ranked-gene-tree-file-name> <gene-tree-topology-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more gene trees with specified branch lengths in the Newick format. The taxon names of gene trees should match the taxon names of the corresponding species tree.   
* ```<gene-tree-topology-file-name>``` (optional input file) contains corresponding gene tree topologies (see Examples). 
* The program outputs a species tree topology in *STtopo.txt* and gene tree topologies along with corresponding probabilities (if input file is given) in *outRankGT.txt*.
  
```
./ranked -uprob <species-tree-file-name> <unranked-gene-tree-file-name>
```
* ```<unranked-gene-tree-file-name>``` contains one unranked gene trees without branch lengths in the Newick format. The taxon names of gene trees should match the taxon names of the corresponding species tree. The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology.   
* The program outputs a species tree topology in *STtopo.txt* and probabilities along with topologies in *outUnrGT.txt*.

## Examples
All files used below can be found in the *tests* folder. 
```
../ranked -rprob st_5taxon.txt rgt_5taxon.txt gtopos_5taxon.txt
```
the output:
```
Total: 0.146615
```
The program also outputs 
```
0.0687959	BE-2-ACD-3-CD-4-
0.0685643	ACD-2-BE-3-CD-4-
0.00925435	ACD-2-CD-3-BE-4-
```

