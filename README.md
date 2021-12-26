# same-class-number

This paper includes code related to the papers:

- "The relative class number one problem for function fields, I" by Kiran S. Kedlaya;
- "The relative class number one problem for function fields, II" by Kiran S. Kedlaya;
- one or more sequel papers in preparation.

Most of the code is packaged in Jupyter notebooks; these are intended to be run using SageMath (tested using version 9.4). When indicated, there is also an external dependency on Magma (tested using version 2.25-5).

In the subdirectory "Shared":
- `cyclic_covers.sage`: Sage subroutines for finding cyclic covers of function fields. Uses Magma.
- `splitting_sequences.sage`: Sage subroutines for computing splitting sequences for finite extensions of function fields.
- `weil_poly_utils.sage`: Sage subroutines for managing Weil polyonmials. This includes converting back and forth between Weil polynomials, point counts, Frobenius traces, and LMFDB labels for isogeny classes of abelian varieties; it also includes calling Sage's Weil polynomials iterator under specific conditions (e.g., with specified Frobenius traces).
- `polys.xlsx`: Excel spreadsheet listing candidate Weil polynomials for purely geometric extensions of function fields over F_2, and indicating which of these correspond to cyclic covers (when known).

In the main subdirectory:

- `orbits.sage`: Sage subroutines for computing orbit representatives for the action of a finite group on subsets of a finite set.

In the subdirectory "Part I":

- `Linear programming.ipynb`: Compute bounds on numbers of rational points on curves over finite fields using the linear programming method.
- `Relative class number 1 for q>2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F_q for q>2, and search for cyclic covers among these. Uses Magma.
- `Weil polynomial bound for q=2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F_2. Updates `polys.xlsx`.
- `Genus 5 hyperelliptic or trigonal.ipynb`: Search for hyperelliptic/trigonal function fields of genus 5 occurring as extensions of genus 1 function fields with relative class number 1. Uses Magma.
- `Genus 5 generic.ipynb`: Search for generic (not hyperelliptic or trigonal) function fields of genus 5 occurring as extensions of genus 1 function fields with relative class number 1. Uses Magma.
- `Cyclic covers for q=2.ipynb`: Search for cyclic covers among candidate Weil polynomials listed in `polys.xlsx`. Uses Magma. Requires and updates `polys.xlsx`.

In the subdirectory "Part II":
- `Non-Galois covers for q=2.ipynb`: Perform various computations needed to rule out the existence of noncyclic covers over F_2. Uses Magma.Requires `polys.xlsx`. 
 
In the subdirectory "Part III":
- To be updated.