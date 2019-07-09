# PRANC
can be used to compute the probabilities of ranked or unranked phylogenetic gene tree topologies given a species tree under the coalescent process.  

## Installation
After downloading the source code, type
```
make
```
This will create an executable called *pranc*, which can be run with some input options listed below.

## Usage
Program options:

| Option        | Input files   | Output files                   |
| ------------- |:-------------| :------------------------------|
| -rprob        | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li><li> file containing gene tree topologies (optional)</li></ul>|<ul><li>STtopo.txt</li><li>outRankGT.txt</li></ul>|
| -uprob        | <ul><li>species tree file</li><li>file containing unranked gene trees (without branch lengths)</li></ul>| <ul><li>STtopo.txt</li><li>outUnrGT.txt</li><li>unrGT.txt</li></ul>|
| -sym        | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li></ul>| <ul><li>STtopo.txt</li><li>out_symbolic.txt</li><li>hist_probs.txt</li></ul>|

All input files should be in the Newick format. All trees are treated as rooted binary trees. We assume an ultrametric species tree (leaves of the tree are all equidistant from the root). The taxon names of gene trees should match the taxon names of the corresponding species tree. User can run *PRANC* with either *-rprob*, *-uprob*, or *-sym* option as shown below.  

```
./pranc -rprob <species-tree-file-name> <ranked-gene-tree-file-name> <gene-tree-topology-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more ranked gene trees in the Newick format.
* ```<gene-tree-topology-file-name>``` (optional input file) contains corresponding gene tree topologies (see Examples). 
* The program outputs a species tree topology (*STtopo.txt*) and probabilities of ranked gene tree topologies (*outRankGT.txt*).
  
```
./pranc -uprob <species-tree-file-name> <unranked-gene-tree-file-name>
```
* ```<unranked-gene-tree-file-name>``` contains one or more unranked gene trees in the Newick format. The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology.   
* The program outputs a species tree topology (*STtopo.txt*), probabilties of unranked gene tree topologies (*unrGT.txt*), and probabilities of ranked topologies (*outUnrGT.txt*).

```
./pranc -sym <species-tree-file-name> <ranked-gene-tree-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more ranked gene trees in the Newick format. The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology.   
* The program outputs a species tree topology (*STtopo.txt*), probabilties of ranked histories (*hist_probs.txt*), and symbolic probabilities (*out_symbolic.txt*).

## Examples
All input files used below can be found in the *tests* folder. 

#### Example 1
```
./pranc -rprob st_5taxon.txt rgt_5taxon.txt gtopos_5taxon.txt
```
output:
```
Total: 0.146615
```
STtopo.txt: 
```
ACD-2-BE-3-CD-4-
```
outRankGT.txt (probabilities and ranked topologies):
```
0.0687959	BE-2-ACD-3-CD-4-
0.0685643	ACD-2-BE-3-CD-4-
0.00925435	ACD-2-CD-3-BE-4-
```

#### Example 2
```
./pranc -rprob st_5taxon.txt rgt_5taxon.txt
```
output:
```
Total: 0.146615
```
STtopo.txt: 
```
ACD-2-BE-3-CD-4-
```
outRankGT.txt (probabilities):
```
0.0687959	
0.0685643	
0.00925435	
```

#### Example 3
```
./pranc -uprob st_5taxon.txt unrgt_5taxon.txt 
```
output:
```
Total: 0.146615
```
STtopo.txt: 
```
ACD-2-BE-3-CD-4-
```
outUnrGT.txt (probabilities and ranked topologies):
```
0.0687959	BE-2-ACD-3-CD-4-
0.0685643	ACD-2-BE-3-CD-4-
0.00925435	ACD-2-CD-3-BE-4-
```
unrGT.txt (unranked tree and probability):
```
((B,E),(A,(C,D)));	0.146615
```

#### Example 4
```
./pranc -sym st_5taxon.txt gt_5taxon.txt 
```
output:
```
Total: 0.0687959
```
hist_probs.txt (ranked histories and probabilities):
```
1234	0.000118525
1233	7.12235e-08
...
1112	0.00373918
1111	0.000909714
```
out_symbolic.txt (first block shows the probability of the ranked history *1234*, second block shows the probability of the ranked history *1233*, etc.)
```
 + (exp(-0*(s1-s2))*1/(1) + exp(-1*(s1-s2))*1/(-1))  * 
(exp(-0*(s2-s3))*1/(1) + exp(-1*(s2-s3))*1/(-1))  * 
(exp(-0*(s3-s4))*1/(1) + exp(-1*(s3-s4))*1/(-1))  * 
2/2

 + (exp(-0*(s1-s2))*1/(1) + exp(-1*(s1-s2))*1/(-1))  * 
(exp(-0*(s2-s3))*1/(2) + exp(-1*(s2-s3))*1/(-1) + exp(-2*(s2-s3))*1/(2))  * 
(exp(-1*(s3-s4))*1/(1))  * 
2/2
...
```
