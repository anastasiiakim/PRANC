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
| :------------------------ |:-----------------------------------------------------------------|:---------------| :-------------------|
| -rprob             |calculates probabilities of ranked gene tree topologies           | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li><li> file containing gene tree topologies (optional)</li></ul>|<ul><li>outRankGT.txt</li></ul>|
| -uprob             |calculates probabilities of unranked gene tree topologies         | <ul><li>species tree file</li><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)</li></ul>| <ul><li>outEachRankTopo.txt</li><li>outUnrGT.txt</li></ul>|
| -sym               |outputs symbolic probabilities of ranked gene tree topologies     | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li></ul>| <ul><li>outSymbolic.txt</li><li>outHistProbs.txt</li></ul>|
| -rtopo             |outputs ranked tree topologies and frequencies of the topologies  | <ul><li>file containing ranked trees (with branch lengths specified)</li></ul>| <ul><li>outRankTopos.txt</li><li>outRankFreqs.txt</li></ul>|
| -utopo             |outputs unranked tree topologies and frequencies of the topologies| <ul><li>file containing unranked trees (without branch lengths; branch lengths will be ignored if given)</li></ul>| <ul><li>outUnrTopos.txt</li><li>outUnrFreqs.txt</li></ul>|
| -write             |outputs tree with ranks instead of branch lengths                 | <ul><li>file containing one ranked gene tree (with branch lengths)</li></ul>|<ul><li>outRankTree.txt</li></ul>|
| -rank_trees        |outputs all ranked topologies that share same unranked topology   | <ul><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)</li></ul>|<ul><li>outRankTopos.txt</li></ul>|
| -mac               |outputs species tree MAC score                                    | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)|<ul><li>outMacScore.txt</li></ul>|
| -cons              |outputs greedy consensus tree without branch lengths              | <ul><li>file containing unranked gene trees (ranked trees will be treated as unranked treees)|<ul><li>outGreedyCons.txt</li></ul>|
| -like_nonni <ul>-rgt</ul>|calculates ML interval lengths of a given species tree topology   | <ul><li>species tree file</li><li>file containing ranked gene trees (with branch lengths)</li></ul>|<ul><li>outNoNniMLTopo.txt</li></ul>|
| -like_nonni <ul>-ugt</ul>|calculates ML interval lengths of a given species tree topology   | <ul><li>species tree file</li><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)|<ul><li>outNoNniMLTopo.txt</li></ul>|
| -like_nni <ul>-rgt</ul>  |estimates ML species tree given a starting tree                   | <ul><li>starting species tree file</li><li>file containing ranked gene trees (with branch lengths)|<ul><li>outWithNniMLTopo.txt</li></ul>|
| -like_nni <ul>-ugt</ul>  |estimates ML species tree given a starting tree                   | <ul><li>starting species tree file</li><li>file containing unranked gene trees (without branch lengths; branch lengths will be ignored if given)|<ul><li>outWithNniMLTopo.txt</li></ul>|



All input files should be in the Newick format. All trees are treated as rooted binary trees. We assume an ultrametric species tree (leaves of the tree are all equidistant from the root). The taxon names of gene trees should match the taxon names of the corresponding species tree. User can run *PRANC* as shown below.  

```
./pranc -rprob <species-tree-file-name> <ranked-gene-tree-file-name> <gene-tree-topology-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more ranked gene trees in the Newick format.
* ```<gene-tree-topology-file-name>``` (optional input file) contains corresponding gene tree topologies (see Examples). 
* The program outputs probabilities of ranked gene tree topologies (*outRankGT.txt*).
  
```
./pranc -uprob <species-tree-file-name> <unranked-gene-tree-file-name>
```
* ```<unranked-gene-tree-file-name>``` contains one or more unranked gene trees in the Newick format.   
* The program ranks an unranked tree and computes probabilities of corresponding ranked gene trees that share the same unranked topology. It outputs probabilties of unranked gene tree topologies (*outUnrGT.txt*) and probabilities of ranked topologies (*outEachRankTopo.txt*).

```
./pranc -sym <species-tree-file-name> <ranked-gene-tree-file-name>
```
* ```<ranked-gene-tree-file-name>``` contains one or more ranked gene trees in the Newick format. 
* The program outputs probabilties of ranked histories (*outHistProbs.txt*) and symbolic probabilities (*outSymbolic.txt*).

```
./pranc -rtopo <ranked-tree-file-name>
```
* ```<ranked-tree-file-name>``` contains one or more ranked trees in the Newick format. 
* The program outputs ranked tree topologies (*outRankTopos.txt*) and frequencies of the topologies (*outRankFreqs.txt*). Species names are separated by "|" character, and groups of clades separated by "-" character followed by clade's rank.
```
./pranc -utopo <unranked-tree-file-name>
```
* ```<unranked-tree-file-name>``` contains one or more unranked trees in the Newick format. 
* The program outputs unranked tree topologies (*outUnrTopos.txt*) and frequencies of the topologies (*outUnrFreqs.txt*). Species names are separated by "|" character, and groups of clades separated by "-" character.


## Examples
All input files used below can be found in the *BIN* folder. 

#### Example 1 (-rprob)
```
./pranc -rprob st_5taxon.txt rgt_5taxon.txt gtopos_5taxon.txt
```
output:
```
Total: 0.146615
```
outRankGT.txt (probabilities and ranked topologies):
```
0.0687959	BE-2-ACD-3-CD-4-
0.0685643	ACD-2-BE-3-CD-4-
0.00925435	ACD-2-CD-3-BE-4-
```

#### Example 2 (-rprob)
```
./pranc -rprob st_5taxon.txt rgt_5taxon.txt
```
output:
```
Total: 0.146615
```
outRankGT.txt (probabilities):
```
0.0687959	
0.0685643	
0.00925435	
```

#### Example 3 (-uprob)
```
./pranc -uprob st_5taxon.txt unrgt_5taxon.txt 
```
output:
```
Total: 0.146615
```
outEachRankTopo.txt (probabilities and ranked topologies):
```
0.0687959	BE-2-ACD-3-CD-4-
0.0685643	ACD-2-BE-3-CD-4-
0.00925435	ACD-2-CD-3-BE-4-
```
outUnrGT.txt (unranked tree and probability):
```
((B,E),(A,(C,D)));	0.146615
```

#### Example 4 (-sym)
```
./pranc -sym st_5taxon.txt gt_5taxon.txt 
```
output:
```
Total: 0.0687959
```
outHistProbs.txt (ranked histories and probabilities):
```
1234	0.000118525
1233	7.12235e-08
...
1112	0.00373918
1111	0.000909714
```
outSymbolic.txt (first block shows the probability of the ranked history *1234*, second block shows the probability of the ranked history *1233*, etc.)
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


#### Example 5 (-rtopo)
```
./pranc -rtopo 5taxa_trees.txt
```

output:

outRankTopos.txt:
```
t1|t2|t3|t4|-2-t1|t3|t4|-3-t1|t4|-4-
t1|t2|t5|-2-t3|t4|-3-t1|t2|-4-
t1|t2|t5|-2-t3|t4|-3-t1|t2|-4-
t2|t3|t4|t5|-2-t2|t5|-3-t3|t4|-4-
t2|t3|t4|t5|-2-t3|t4|-3-t2|t5|-4-
```
outRankFreqs.txt: 
```
2	t1|t2|t5|-2-t3|t4|-3-t1|t2|-4-
1	t2|t3|t4|t5|-2-t3|t4|-3-t2|t5|-4-
1	t2|t3|t4|t5|-2-t2|t5|-3-t3|t4|-4-
1	t1|t2|t3|t4|-2-t1|t3|t4|-3-t1|t4|-4-
```


#### Example 6 (-utopo)
```
./pranc -utopo 5taxa_trees.txt
```
output:

outUnrTopos.txt:
```
t1|t4|-t1|t3|t4|-t1|t2|t3|t4|-t1|t2|t3|t4|t5|-
t1|t2|-t3|t4|-t1|t2|t5|-t1|t2|t3|t4|t5|-
t1|t2|-t3|t4|-t1|t2|t5|-t1|t2|t3|t4|t5|-
t2|t5|-t3|t4|-t2|t3|t4|t5|-t1|t2|t3|t4|t5|-
t2|t5|-t3|t4|-t2|t3|t4|t5|-t1|t2|t3|t4|t5|-
```
outUnrFreqs.txt: 
```
2	t2|t5|-t3|t4|-t2|t3|t4|t5|-t1|t2|t3|t4|t5|-
2	t1|t2|-t3|t4|-t1|t2|t5|-t1|t2|t3|t4|t5|-
1	t1|t4|-t1|t3|t4|-t1|t2|t3|t4|-t1|t2|t3|t4|t5|-
```


#### Example 7 (-utopo)
```
./pranc -utopo unrgts.txt
```

output:

outUnrTopos.txt:
```
t3|t4|-t6|t7|-t1|t6|t7|-t3|t4|t5|-t1|t6|t7|t8|-t1|t2|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t7|t8|-t1|t2|t3|-t6|t7|t8|-t1|t2|t3|t4|-t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t1|t2|-t3|t4|-t5|t6|-t7|t8|-t5|t6|t7|t8|-t3|t4|t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
t3|t4|-t6|t7|-t1|t6|t7|-t3|t4|t5|-t1|t6|t7|t8|-t1|t2|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
```
outUnrFreqs: 
```
3	t1|t2|-t3|t4|-t7|t8|-t1|t2|t6|-t5|t7|t8|-t3|t4|t5|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
2	t3|t4|-t6|t7|-t1|t6|t7|-t3|t4|t5|-t1|t6|t7|t8|-t1|t2|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
1	t1|t2|-t7|t8|-t1|t2|t3|-t6|t7|t8|-t1|t2|t3|t4|-t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
1	t1|t2|-t3|t4|-t5|t6|-t7|t8|-t5|t6|t7|t8|-t3|t4|t5|t6|t7|t8|-t1|t2|t3|t4|t5|t6|t7|t8|-
```

#### Example 8 (-write)
```
./pranc -write st_5taxon.txt 
```
output:

outRankTree.txt:
```
((B:2,E:2):2,(A:3,(C:1,D:1):2):1);
```

#### Example 9 (-rank_trees)
```
./pranc -rank_trees unrgt_5taxon.txt
```
output:

outRankTopos.txt:
```
((B:3,E:3):1,(A:2,(C:1,D:1):1):2);
((B:2,E:2):2,(A:3,(C:1,D:1):2):1);
((B:1,E:1):3,(A:3,(C:2,D:2):1):1);
```

#### Example 10 (-mac)
```
./pranc -mac st_5taxon.txt rgt_5taxon.txt
```
output:

outMacScore.txt:
```
2
```

#### Example 11 (-cons)
```
./pranc -cons unrgts.txt
```
output:

outGreedyCons.txt:
```
((t6,(t1,t2)),((t3,t4),(t5,(t7,t8))));
```

#### Example 12 (-like_nonni -rgt)
```
./pranc -like_nonni st_5taxon.txt -rgt rgt_5taxon.txt
```
output:

outNoNniMLTopo.txt (your estimated branch lengths will be slightly different):
```
((B:0.375797,E:0.375797):7.173825,(A:1.549622,(C:0.100000,D:0.100000):1.449622):6.000000);
```

#### Example 12 (-like_nonni -ugt)
```
./pranc -like_nonni st_5taxon.txt -ugt ugt_5taxon.txt
```
output:

outNoNniMLTopo.txt (your estimated branch lengths will be slightly different):
```
((B:0.100100,E:0.100100):0.000200,(A:0.100200,(C:0.100000,D:0.100000):0.000200):0.000100);
```

#### Example 13 (-like_nni -rgt)
```
./pranc -like_nni st_5taxon.txt -rgt rgt_5taxon.txt
```
output:

outWithNniMLTopo.txt (your estimated branch lengths and topology might be slightly different):
```
((B:0.376472,E:0.376472):7.173572,(A:1.550044,(C:0.100000,D:0.100000):1.450044):6.000000);
```

#### Example 14 (-like_nni -ugt)
```
./pranc -like_nni st_5taxon.txt -ugt ugt_5taxon.txt
```
output:

outWithNniMLTopo.txt (your estimated branch lengths and topology might be slightly different):
```
(((E:0.100000,D:0.100000):0.592201,A:0.692201):0.255117,(B:0.246954,C:0.246954):0.700365);
```
