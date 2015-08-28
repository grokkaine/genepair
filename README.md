# genepair
Measuring gene co-expression by means of mutual information

Read the notes.txt for usage instructions.

TODO:
Acknowledge Carrera et al.

In the (arguably rare) case that the number of input genes is larger than 46341, the `calloc()` call in `computeCLR()` will fail as the integer overflows. A simple fix is to change the type of ngenes from `int` to `size_t`.

gcc -L/usr/bin -lz -lm genepair.c genepair_f.c -o genepair
On Ubunt, lm must be placed after the code:
gcc -L/usr/bin -lz genepair.c genepair_f.c -lm -o genepair

Mutual information (mi) computation:
nohup ./genepair 6 microarray_file.txt mi_file.txt ngenes nexp fromgene togene num_bins spline_order &

Parameters:
6 - indicates that MI values must be computed (use 1 when there are no undefined values, it's a little faster)
microarray_file.txt - Must be space separated, the number of lines must be equal to the number of genes, each line has a list of space separated values equal to the number of experiments.
(the size of the microarray matrix is ngenes*nexp)
mi_file.txt - The program will write here the MI values, it must be a lower triangular matrix of gene-correlaiton pairs, with the diagonal excluded (thus the number of lines must be smaller by 1 then the number of genes)
(values are appended to this file, so when restarting the computation, the file must be deleted or renamed if it already exists)
(after the file is computed it must be checked if the number of lines is correct)
ngenes - The number of genes in the microarray
nexp - The number of experiments in the microarray
fromgene - 0
togene - ngenes
numbins - (7) The number of bins used for computing the probability distribution of gene expression. Normally it would not be related to the number of experiments and the expression value range, but we use a constant value of 7 in practice.
spline_order - (3) The order of the splines used to smooth the probability distribution functions (we use 3 by default).

When it is over and the mi_file is generated do next step:
CLR computation (Context Likelihood of Relatedness - This is basically a way to measure the local significance of every gene correlation)
nohup ./genepair 2 mi_file.txt ngenes clr_file.txt inputfiletype &
(takes much less time and memory to run than the MI computation)

Parameters:
2 - indicates that CLR must be computed
mi_file.txt - Input file of MI values, must be in lower triangular form, and must be properly checked before running the program (see above). As long as the value are space separated, the program will not check for errors in the imput file, it is you who must ensure that!
ngenes - The number of genes in the microarray
clr_file.txt - The outputted CLR file, must always be checked if the number of lines is correct. The format is the same as mi_file.txt, a lower matrix with number of lines must be smaller by 1 then the number of genes.
inputfiletype - (1) Must be 1 to indicate that the type if the mi_file is a lower triangular matrix.

To get the number of experiments:
sed -n 1p name_microarray.txt |wc -w

Get the number of genes:
wc -l name_microarray.txt
