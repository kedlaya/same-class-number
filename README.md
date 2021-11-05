# same-class-number
This paper includes code related to the papers "The relative class number one problem for function fields, I" by Kiran S. Kedlaya and "The relative class number one problem for function fields, II" (in progress). Most of the code is packaged in Jupyter notebooks; these are intended to be run using SageMath (tested using version 9.5beta2) and Magma (tested using version 2.25-5).

- `auxiliary.sage`: Shared subroutines used by multiple notebooks.
- `Linear programming.ipynb`: Compute bounds on numbers of rational points on curves over finite fields using the linear programming method.
- `Relative class number 1 for q>2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F_q for q>2, and search for cyclic covers among these.
- `Weil polynomial bound for q=2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F_2. Updates `polys.xlsx`.
- `Cyclic covers for q=2.ipynb`: Search for cyclic covers among candidate Weil polynomials listed in `polys.xlsx`. Updates `polys.xlsx`.
- `polys.xlsx`: Output file listing candidate Weil polynomials for purely geometric extensions of function fields over F_2, and indicating which of these correspond to cyclic covers.
- `Non-Galois covers for q=2.ipynb`: Perform various computations needed to rule out the existence of noncyclic covers over F_2.
