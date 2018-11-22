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

## Examples

