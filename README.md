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
$$
\frac{n!}{k!(n-k)!} = {n \choose k}
$$
```

```
\begin{tikzpicture}[scale=0.6, transform shape, every label/.append style={font=\LARGE},
         rec/.style={shape=circle, draw=red,fill=white, line width=1},
         regr/.style={shape=circle, fill=black, draw=black, inner sep=1pt,minimum size=1pt},
         lsq/.style={thick,->}]

     \draw(1,1)node[label=below:{A}](x1){};
     \draw(3,1)node[label=below:{B}](x2){};
     \draw(4.5,1)node[label=below:{C}](x3){};
     \draw(7,1)node[label=below:{D}](x4){};
     \draw(8.5,1)node[label=below:{E}](x5){};


     \draw[-,black](0.5,1)  -- (1.5,1); %bottomA
     \draw[-,black](0.5,1)  -- (1,3) -- (1.25,4) -- (1.4,4.5) -- (2.4,6.272728) -- (2.4,7.272728); %leftA

     \draw[-,black](1.5,1)  -- (2,3); %rightA
     \draw[-,black](2,3)  -- (2.5,1); %leftB
     \draw[-,black](2.5,1)  -- (3.5,1); %bottomB
     \draw[-,black](3.5,1)  -- (3,3); %bottom part of rightB


     \draw[-,black](4,1)  -- (5,1); %bottomC
     \draw[-,black](4,1)  -- (3,4.5); %leftC
     \draw[-,black](3,4.5)  -- (2.75,4) -- (3,3); %upper part of rightB
     \draw[-,black](5,1)  -- (4,4.5) -- (3.63636364,5.77272728) -- (4.38636364,6.272728); %rightC
     \draw[-,black](4,4.5); %rightC

     \draw[-,black](6.875,2.5) -- (6.5,1)  -- (7.5,1) -- (7.75,2); %pants D
     \draw[-,black](7.75,2) -- (8,1)  -- (9,1) -- (8.625,2.5); %pants E

    \draw[-,black](6.875,2.5) -- (4.38636364,6.272728); %leftD
    \draw[-,black](8.625,2.5) -- (6.13636364,6.272728) -- (6.13636364,7.272728); %rightE



    \draw(0.5,2)node[label=left:{$s_4$}](x6){};
    \draw(0.5,3)node[label=left:{$s_3$}](x7){};
    \draw(0.5,4.5)node[label=left:{$s_2$}](x8){};
    \draw(0.5,6.272728)node[label=left:{$s_1$}](x9){};

    \draw[dotted,blue,thick](0.5,2)  -- (9.5,2);
    \draw[dotted,blue,thick](0.5,3)  -- (9.5,3);
    \draw[dotted,blue,thick](0.5,4.5)  -- (9.5,4.5);
    \draw[dotted,blue,thick](0.5,6.272728)  -- (9.5,6.272728);

    \draw(9,2.5)node[label=right:{$\tau_4$}](x10){};
    \draw(9,3.75)node[label=right:{$\tau_3$}](x11){};
    \draw(9,5.386364)node[label=right:{$\tau_2$}](x11){};
    \draw(9,6.772728)node[label=right:{$\tau_1$}](x11){};
    
  
    \draw(6.4,4.8)node[regr,label=left:{$u_{4}$}](u23){};
    \draw(2.5,5.2)node[regr,label=right:{$u_{3}$}](u22){};
    \draw(3.2,6)node[regr,label=left:{$u_{2}$}](u21){};
    \draw(4.4,7)node[regr,label=above:{$u_{1}$}](u11){};


    \draw(6.8,1)  to [bend right = 10](7,1.5) to [bend left = 10](7.1,2) to [bend right = 10](7.1,2.5) to [bend left = 10](6.8,3) to [bend right = 10](6.6,3.5) to [bend left = 10](6.3,4)  to [bend right = 10](6.25,4.4) to [bend left = 10](u23); %bottom left line of u23
    \draw(8.5,1)  to [bend right = 10](8.3,1.5) to [bend left = 10](8.2,2) to [bend right = 10](8.1,2.5) to [bend left = 10](7.9,3) to [bend right = 10](7.7,3.5) to [bend left = 10](7.3,4) to [bend right = 10](7,4.4) to  [bend left = 10](u23); %bottom right line of u23
    \draw(u23)  to [bend right = 10](6,5.3) to [bend left = 10](5.7,5.8) to [bend right = 10](5.2,6.5) to [bend left = 10](4.9,6.75) to [bend right = 10](u11); %upper line of u23


    \draw(1,1)  to [bend right = 10](1.2,1.5) to [bend left = 10](1.4,2) to [bend right = 10](1.6,2.5) to [bend left = 10](1.7,3) to [bend right = 10](1.8,3.5) to [bend left = 10](2,4)  to [bend right = 10](2.1,4.4) to [bend left = 10](u22); %bottom left line of u22

    \draw(3.2,1)  to [bend right = 10](3.1,1.5) to [bend left = 10](2.95,2) to [bend right = 10](2.9,2.5) to [bend left = 10](2.8,3) to [bend right = 10](2.7,3.5) to [bend left = 10](2.55,4) to [bend right = 10](2.6,4.4) to  [bend left = 10](u22); %bottom right line of u22

    \draw(u22)  to [bend right = 10](2.7,5.4) to [bend left = 10](u21) to [bend right = 10](3.7,6.4) to [bend left = 10](4.1,6.75) to [bend right = 10](u11); %upper line of u22

    \draw(4.55,1)  to [bend right = 10](4.5,1.5) to [bend left = 10](4.3,2) to [bend right = 10](4.2,2.5) to [bend left = 10](4.1,3) to [bend right = 10](4,3.5) to [bend left = 10](3.9,4) to [bend right = 10](3.8,4.3) to [bend left = 10](3.5,5) to [bend right = 10](u21); %bottom line of u21
  \node at (5.5, 11) {\huge \textbf{(a)} Matching ranked gene tree};
   \node at (5.8, 10) {\huge Ranked history (1,2,2,2)};
       \end{tikzpicture}

```


## Examples

