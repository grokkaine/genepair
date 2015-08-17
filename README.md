# genepair
Measuring gene co-expression by means of mutual information

TODO:
Acknowledge Carrera et al.

In the (arguably rare) case that the number of input genes is larger than 46341, the `calloc()` call in `computeCLR()` will fail as the integer overflows. A simple fix is to change the type of ngenes from `int` to `size_t`.
