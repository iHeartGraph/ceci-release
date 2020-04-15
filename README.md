# ceci-release
Source Code for "CECI: Compact Embedding Cluster Index for Scalable Subgraph Matching"

## Introduction
CECI is a subgraph matching system that works for undirected query/data graphs with vertex labels. It uses subgraph isomorphism as the embedding syntax and utilizes multiple threads to speed up the listing process. 

## Data and Query Graph Format
The graphs are represented with vertex and edge lists concatenated in single file. Following example shows file representing a Triangle i.e., qg1.

```markdown
t # 0
v 0 -1
v 1 -1
v 2 -1
e 0 1 0
e 1 2 0
e 2 0 0
```
The line starting with 't' represents the graph identifier. For a file containing single graph, this can simply be used as it is.

The line starting with 'v' are the vertices of the graph. The first number following 'v' is vertex identifier and the second number is identifier for vertex label. 

The line starting with 'e' are the edges of the graph. The remaining three numbers are first endpoint, second endpoint and label of the edge respectively. The edges are undirected, which means graph traversal can occur in bothe direction between two endpoints. An expmple data graph (dg) is shown below.
```markdown
t # 0
v 0 0
v 1 0
v 2 0
v 3 0
v 4 0
e 0 1 0
e 0 2 0
e 0 3 0
e 0 4 0
e 1 2 0
e 1 3 0
e 1 4 0
e 2 3 0
e 2 4 0
e 3 4 0
```

## Running CECI 
Once you compile the CECI, it should result in a binary named **ceci**. **ceci** only needs 2 additional inputs i.e. **a data graph** and **a query graph**.
```markdown
./ceci ./dg ./qg1 
```
Sample datasets, queries and scripts are included in the repository.

For any help, enquiries and comments, please contact me at bhattarai_b@gwu.edu 
