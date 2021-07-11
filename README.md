# wavemodes

Lamb mode eigenvalue computation following method of

> Pagneux, V., & Maurel, A. (2001). Determination of Lamb mode eigenvalues. The Journal of the Acoustical Society of America, 110(3), 1307-1314. https://doi.org/10.1121/1.1391248

A PDF of the article can be downloaded freely on [ResearchGate](https://www.researchgate.net/publication/11776614_Determination_of_Lamb_mode_eigenvalues).

The script `copper.m` located in the current folder attempts to reproduce Fig. 4a of the original article.

## Basis functions

[Symmetric basis functions](figs/sym.png)
[Antisymmetric basis functions](figs/asym.png)

The plot of the basis functions was produced using

```gnuplot
n = 5
h = 0.01
a(n) = (n-1)*pi/h
b(n) = (n-0.5)*pi/h

set terminal pngcairo enhanced

set output "sym.png"
plot [-h:h] for [i=1:n] cos(a(i)*x) lw 3 t sprintf('n = %d',i)

set output "ssym.png"
plot [-h:h] for [i=1:n] sin(b(i)*x) lw 3 t sprintf('n = %d',i)
```

## Fortran code

A Fortran implementation of the wave mode calculator is in `lamb.f90`. To compile the code use:

```
gfortran -O3 lamb.f90 -o lamb -llapack -lblas
```

Running the code with the command `$ lamb` will produce a file `eigvals.txt` that can be plotted using gnuplot.