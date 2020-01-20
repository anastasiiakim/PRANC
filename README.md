# PRANC
can be used to compute the probabilities of ranked or unranked phylogenetic gene tree topologies given a species tree under the coalescent process. It can also output "democratic vote" (most frequent) ranked or unranked topologies. *PRANC* can also estimate the maximum likelihood species tree with branch lengths from the sample of ranked or unranked gene tree topologies. Greedy consensus tree can be used as a starting tree. Also, trees selected by minimization of ancient coalescence (MAC) criterion can be used as starting trees.

## Installation
After downloading the source code, go to SRC directory and type
```
make
```
This will create an executable called *pranc*, which can be run from BIN with some input options listed below.

## Usage
Program options:

| Option             |Description                                                       | Input files    | Output files        |
| ------------------ |:-----------------------------------------------------------------|:---------------| :-------------------|
| -rprob             |calculates probabilities of ranked gene tree topologies           | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li><li> file containing gene tree topologies (optional)</li></ul>|<ul><li>outRankGT.txt</li></ul>|
| -uprob             |calculates probabilities of unranked gene tree topologies         | <ul><li>species tree file</li><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)</li></ul>| <ul><li>outEachRankTopo.txt</li><li>outUnrGT.txt</li></ul>|
| -sym               |outputs symbolic probabilities of ranked gene tree topologies     | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li></ul>| <ul><li>outSymbolic.txt</li><li>outHistProbs.txt</li></ul>|
| -rtopo             |outputs ranked tree topologies and frequencies of the topologies  | <ul><li>file containing ranked trees (with branch lengths specified)</li></ul>| <ul><li>outRankTopos.txt</li><li>outRankFreqs.txt</li></ul>|
| -utopo             |outputs unranked tree topologies and frequencies of the topologies| <ul><li>file containing unranked trees (without branch lengths; branch lengths will be ignored if given)</li></ul>| <ul><li>outUnrTopos.txt</li><li>outUnrFreqs.txt</li></ul>|
| -write             |outputs tree with ranks instead of branch lengths                 | <ul><li>file containing one ranked gene tree (with branch lengths)</li></ul>|<ul><li>outRankTree.txt</li></ul>|
| -rank_trees        |outputs all ranked topologies that share same unranked topology   | <ul><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)</li></ul>|<ul><li>outRankTopos.txt</li></ul>|
| -acmin             |outputs species tree MAC score                                    | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)|<ul><li>outMacScore.txt</li></ul>|
| -cons              |outputs greedy consensus tree without branch lengths              | <ul><li>file containing unranked gene trees (ranked trees will be treated as unranked treees)|<ul><li>outGreedyCons.txt</li></ul>|
| -like_nonni -rgt   |calculates ML interval lengths of a given species tree topology   | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li></ul>|<ul><li>outNoNniMLTopo.txt</li></ul>|
| -like_nonni -ugt   |calculates ML interval lengths of a given species tree topology   | <ul><li>species tree file</li><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)|<ul><li>outNoNniMLTopo.txt</li></ul>|
| -like_nni -rgt     |estimates ML species tree given a starting tree                   | <ul><li>starting species tree file</li><li>file containing ranked gene trees (with branch lengths)|<ul><li>outWithNniMLTopo.txt</li></ul>|
| -like_nni -ugt     |estimates ML species tree given a starting tree                   | <ul><li>starting species tree file</li><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)|<ul><li>outWithNniMLTopo.txt</li></ul>|



All input files should be in the Newick format. All trees are treated as rooted binary trees. We assume an ultrametric species tree (leaves of the tree are all equidistant from the root). The taxon names of gene trees should match the taxon names of the corresponding species tree. User can run *PRANC* with either *-rprob*, *-uprob*, *-sym*, *-rtopo*, or *-utopo* option as shown below.  

```
./pranc -rprob <species-tree-file-name> <ranked-gene-tree-file-name> <gene-tree-topology-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more ranked gene trees in the Newick format.
* ```<gene-tree-topology-file-name>``` (optional input file) contains corresponding gene tree topologies (see Examples). 
* The program outputs a species tree topology (*STtopo.txt*) and probabilities of ranked gene tree topologies (*outRankGT.txt*).
  
```
./pranc -uprob <species-tree-file-name> <unranked-gene-tree-file-name>
```
* ```<unranked-gene-tree-file-name>``` contains one or more unranked gene trees in the Newick format.   
* The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology. It outputs a species tree topology (*STtopo.txt*), probabilties of unranked gene tree topologies (*unrGT.txt*), and probabilities of ranked topologies (*outUnrGT.txt*).

```
./pranc -sym <species-tree-file-name> <ranked-gene-tree-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more ranked gene trees in the Newick format. 
* The program outputs a species tree topology (*STtopo.txt*), probabilties of ranked histories (*hist_probs.txt*), and symbolic probabilities (*out_symbolic.txt*).

```
./pranc -rtopo <ranked-tree-file-name>
```
* ```<ranked-tree-file-name>``` contains one or more ranked trees in the Newick format. 
* The program outputs ranked tree topologies (*rtopos.txt*) and frequencies of the topologies (*rfreqs.txt*). Species names are separated by "|" character, and groups of clades separated by "-" character followed by clade's rank.
```
./pranc -utopo <unranked-tree-file-name>
```
* ```<unranked-tree-file-name>``` contains one or more unranked trees in the Newick format. 
* The program outputs unranked tree topologies (*utopos.txt*) and frequencies of the topologies (*ufreqs.txt*). Species names are separated by "|" character, and groups of clades separated by "-" character.


## Examples
All input files used below can be found in the *BIN* folder. 

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


#### Example 5
```
./pranc -rtopo 5taxa_trees.txt
```

output:

rtopos.txt:
```
t1|t2|t3|t4|-2-t1|t3|t4|-3-t1|t4|-4-
t1|t2|t5|-2-t3|t4|-3-t1|t2|-4-
t1|t2|t5|-2-t3|t4|-3-t1|t2|-4-
t2|t3|t4|t5|-2-t2|t5|-3-t3|t4|-4-
t2|t3|t4|t5|-2-t3|t4|-3-t2|t5|-4-
```
rfreqs.txt: 
```
2	t1|t2|t5|-2-t3|t4|-3-t1|t2|-4-
1	t2|t3|t4|t5|-2-t3|t4|-3-t2|t5|-4-
1	t2|t3|t4|t5|-2-t2|t5|-3-t3|t4|-4-
1	t1|t2|t3|t4|-2-t1|t3|t4|-3-t1|t4|-4-
```


#### Example 6
```
./pranc -utopo 5taxa_trees.txt
```
output:

utopos.txt:
```
t1|t4|-t1|t3|t4|-t1|t2|t3|t4|-t1|t2|t3|t4|t5|-
t1|t2|-t3|t4|-t1|t2|t5|-t1|t2|t3|t4|t5|-
t1|t2|-t3|t4|-t1|t2|t5|-t1|t2|t3|t4|t5|-
t2|t5|-t3|t4|-t2|t3|t4|t5|-t1|t2|t3|t4|t5|-
t2|t5|-t3|t4|-t2|t3|t4|t5|-t1|t2|t3|t4|t5|-
```
ufreqs.txt: 
```
2	t2|t5|-t3|t4|-t2|t3|t4|t5|-t1|t2|t3|t4|t5|-
2	t1|t2|-t3|t4|-t1|t2|t5|-t1|t2|t3|t4|t5|-
1	t1|t4|-t1|t3|t4|-t1|t2|t3|t4|-t1|t2|t3|t4|t5|-
```


#### Example 7
```
./pranc -utopo unrgts.txt
```

output:

utopos.txt:
```
t3|t4|-t6|t7|-t1|t6|t7|-t3|t4|t5|-t1|t6|t7|t8|-t1|t2|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t7|t8|-t1|t2|t3|-t6|t7|t8|-t1|t2|t3|t4|-t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t5|t6|-t7|t8|-t5|t6|t7|t8|-t3|t4|t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t3|t4|-t6|t7|-t1|t6|t7|-t3|t4|t5|-t1|t6|t7|t8|-t1|t2|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
```
ufreqs.txt: 
```
3	t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
2	t3|t4|-t6|t7|-t1|t6|t7|-t3|t4|t5|-t1|t6|t7|t8|-t1|t2|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
1	t1|t2|-t7|t8|-t1|t2|t3|-t6|t7|t8|-t1|t2|t3|t4|-t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
1	t1|t2|-t3|t4|-t5|t6|-t7|t8|-t5|t6|t7|t8|-t3|t4|t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
```
