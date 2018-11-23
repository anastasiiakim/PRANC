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
Note: all input files should be in the Newick format. All trees are treates as rooted binary trees. The species tree assumed to be an ultrametric.  

```
<p align="center"><img src="https://rawgit.com/anastasiiakim/ranked_probs/None/svgs/32737e0a8d5a4cf32ba3ab1b74902ab7.svg?invert_in_darkmode" align=middle width=127.89183pt height=39.30498pt/></p>
```
```
<p align="center"><img src="https://rawgit.com/anastasiiakim/ranked_probs/None/svgs/82d10f5787d4ca7bf99a377a7a4c8ef4.svg?invert_in_darkmode" align=middle width=697.94505pt height=35.19318pt/></p>
```

## Examples

