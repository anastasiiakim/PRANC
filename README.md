# PRANC
can be used to compute the probabilities of ranked or unranked phylogenetic gene tree topologies given a species tree under the coalescent process.  

## Installation
After downloading the source code, type
```
make
```
This will create an executable called pranc, which can be futher run with some input options listed below.

## Usage
Program options:

| Option        | Input files   | Output files                   |
| ------------- |:-------------| :------------------------------|
| -rprob        | <ul><li>input a species tree file</li><li>input a file containing ranked gene trees (with branch lengths)</li><li> input a file containing gene tree topologies (optional)</li></ul>|<ul><li>STtopo.txt</li><li>outRankGT.txt</li></ul>|
| -uprob        | <ul><li>input a species tree file</li><li>input a file containing unranked gene trees (without branch lengths)</li></ul>| <ul><li>STtopo.txt</li><li>outUnrGT.txt</li><li>unrGT.txt</li></ul>|
| -sym        | <ul><li>input a species tree file</li><li>input a file containing ranked gene trees (with branch lengths)</li></ul>| <ul><li>STtopo.txt</li><li>out_symbolic.txt</li><li>hist_probs.txt</li></ul>|

All input files should be in the Newick format. All trees are treated as rooted binary trees. We assume an ultrametric species tree (leaves of the tree are all equidistant from the root). User can run *PRANC* with either *-rprob*, *-uprob*, or *-sym* option as shown below.  

```
./pranc -rprob <species-tree-file-name> <ranked-gene-tree-file-name> <gene-tree-topology-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more gene trees with specified branch lengths in the Newick format. The taxon names of gene trees should match the taxon names of the corresponding species tree.   
* ```<gene-tree-topology-file-name>``` (optional input file) contains corresponding gene tree topologies (see Examples). 
* The program outputs a species tree topology in *STtopo.txt* and gene tree topologies (if input file is given) along with corresponding probabilities in *outRankGT.txt*.
  
```
./pranc -uprob <species-tree-file-name> <unranked-gene-tree-file-name>
```
* ```<unranked-gene-tree-file-name>``` contains one or more gene trees with specified branch lengths in the Newick format. The taxon names of gene trees should match the taxon names of the corresponding species tree. The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology.   
* The program outputs a species tree topology in *STtopo.txt*, probabilties of unranked gene tree topologies in *unrGT.txt*, and probabilities along with ranked topologies in *outUnrGT.txt*.

```
./pranc -sym <species-tree-file-name> <ranked-gene-tree-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more ranked gene trees with branch lengths in the Newick format. The taxon names of gene trees should match the taxon names of the corresponding species tree. The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology.   
* The program outputs a species tree topology in *STtopo.txt*, probabilties of ranked histories in *hist_probs.txt*, and symbolic probabilities in *out_symbolic.txt*.

## Examples
All files used below can be found in the *tests* folder. 
```
../pranc -rprob st_5taxon.txt rgt_5taxon.txt gtopos_5taxon.txt
../pranc -uprob st_5taxon.txt unrgt_5taxon.txt 
```
Both options will give the following output:
```
Total: 0.146615
```
STtopo.txt: 
```
ACD-2-BE-3-CD-4-
```
outRankGT.txt/outUnrGT.txt:
```
0.0687959	BE-2-ACD-3-CD-4-
0.0685643	ACD-2-BE-3-CD-4-
0.00925435	ACD-2-CD-3-BE-4-
```
unrGT.txt
```
((B,E),(A,(C,D)));	0.146615
```

```
../pranc -sym st_5taxon.txt gt_5taxon.txt 
```

```
Total: 0.0687959
```

```
1234	0.000118525
1233	7.12235e-08
...
1112	0.00373918
1111	0.000909714
```

```
 + (exp(-0*(s1-s2))*1/(1) + exp(-1*(s1-s2))*1/(-1))  * 
(exp(-0*(s2-s3))*1/(1) + exp(-1*(s2-s3))*1/(-1))  * 
(exp(-0*(s3-s4))*1/(1) + exp(-1*(s3-s4))*1/(-1))  * 
2/2

 + (exp(-0*(s1-s2))*1/(1) + exp(-1*(s1-s2))*1/(-1))  * 
(exp(-0*(s2-s3))*1/(2) + exp(-1*(s2-s3))*1/(-1) + exp(-2*(s2-s3))*1/(2))  * 
(exp(-1*(s3-s4))*1/(1))  * 
2/2
```
