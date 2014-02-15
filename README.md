phylo
=====

Software with few features for calculating trees of maximum parsimony, given an input matrix.

Supports unordered characters (including unknown), outputs tree in newick format.

Search algorithms

- Branch and bound
- Tree bisection and reconnection (TBR)
- Parsimony ratchet

Tree generation

- Random tree generation (dumb.cpp)
- Random stepwise addition (smart.cpp). Slow but generates good starting trees.
 
