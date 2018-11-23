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

The input may look like this:
```
./ranked -rprob <species-tree-file-name> <ranked-gene-tree-file-name>
./ranked -uprob <species-tree-file-name> <unranked-gene-tree-file-name>
```
Note: all input files should be in the Newick format. All trees are treates as rooted binary trees. The species tree assumed to be an ultrametric (leaves of the tree are all equidistant from the root).  

```
<p align="center"><img src="https://github.com/anastasiiakim/ranked-coal/svgs/32737e0a8d5a4cf32ba3ab1b74902ab7.svg?invert_in_darkmode" align=middle width=127.89183pt height=39.30498pt/></p>
```

```
<p align="center"><img src="https://rawgit.com/anastasiiakim/ranked_probs/master/svgs/ed8e99f2b3f7a873819ad793ebc15ae1.svg" align=middle width=421.4991pt height=33.965745pt/></p>

```


## Examples

