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
\begin{tikzpicture}
\newcounter{density}
\setcounter{density}{20}
    \def\couleur{red}
    \path[coordinate] (0,0)  coordinate(A)
                ++( 60:6cm) coordinate(B)
                ++(-60:6cm) coordinate(C);
    \draw[fill=\couleur!\thedensity] (A) -- (B) -- (C) -- cycle;
    \foreach \x in {1,...,15}{%
        \pgfmathsetcounter{density}{\thedensity+10}
        \setcounter{density}{\thedensity}
        \path[coordinate] coordinate(X) at (A){};
        \path[coordinate] (A) -- (B) coordinate[pos=.15](A)
                            -- (C) coordinate[pos=.15](B)
                            -- (X) coordinate[pos=.15](C);
        \draw[fill=\couleur!\thedensity] (A)--(B)--(C)--cycle;
    }
\end{tikzpicture}
```


## Examples

