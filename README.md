# same-class-number

This paper includes code related to the papers "The relative class number one problem for function fields, I, II, III" by Kiran S. Kedlaya. Most of the code is packaged in Jupyter notebooks; these are intended to be run using SageMath (tested using version 9.5beta9). When indicated, there is also an external dependency on Magma (tested using version 2.26-9).

In the subdirectory "Shared":

- `cyclic_covers.sage`: Sage subroutines for finding cyclic covers of function fields. Uses Magma.
- `weil_poly_utils.sage`: Sage subroutines for managing Weil polynomials. This includes converting back and forth between Weil polynomials, point counts, Frobenius traces, and LMFDB labels for isogeny classes of abelian varieties; it also includes calling Sage's Weil polynomials iterator under specific conditions (e.g., with specified Frobenius traces).
- `polys.xlsx`: Excel spreadsheet listing candidate Weil polynomials for purely geometric extensions of function fields over F_2, and indicating which of these correspond to cyclic covers (when known). This is generated and modified by various notebooks.

In the subdirectory "Part I":

- `Linear programming.ipynb`: Compute bounds on numbers of rational points on curves over finite fields using the linear programming method.
- `Relative class number 1 for q>2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F_q for q>2, and search for cyclic covers among these. Uses Magma.
- `Genus 5 data.ipynb`: Process data files from Dragutinović (see above) to build a table of genus 5 curves over F_2 sorted by their Weil polynomials. Creates a (large) temporary file `genus5byweilpoly.txt`.
- `Weil polynomial bound for q=2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F_2. Requires the temporary file `genus5byweilpolyt.txt`. Updates `polys.xlsx`.
- `Cyclic covers for q=2.ipynb`: Search for cyclic covers among candidate Weil polynomials listed in `polys.xlsx`. Uses Magma and data from Dragutinović (https://github.com/DusanDragutinovic/MT_Curves). Requires the temporary file `genus5byweilpolyt.txt`. Requires and updates `polys.xlsx`.
- `HyperellipticCurvesData.txt`, `TrigonalCurvesWithAutomorphisms.txt`, `Complete_Intersections.txt`, `Pts_Count_Complete_Intersections.txt`: data files from Dragutinović (see below).

In the subdirectory "Part II":

- `splitting_sequences.sage`: Sage subroutines for computing splitting sequences for finite extensions of function fields.
- `Non-Galois covers for q=2.ipynb`: Perform various computations needed to rule out the existence of noncyclic covers over F_2. Uses Magma. Requires `polys.xlsx`. 
 
In the subdirectory "Part III":

- `orbits.sage`: Sage subroutines for computing orbit representatives for the action of a finite group on subsets of a finite set.
- `linalg.sage`: Sage subroutines for linear algebra.
- `preamble.sage`: Common declarations for the notebooks in this folder. This includes some initialization code   plus a function `closeout` run at the end to report results.
- The other files in this folder are Jupyter notebooks covering one (or more) Brill-Noether strata in genus 6 or 7. These require `polys.xlsx` and depend on Magma.
  - Generic (non-tetragonal) curves of genus 7 are handled in a series of three notebooks:
     - `Genus 7 generic, part 1.ipynb` performs the orbit lookup tree calculation. Creates a temporary file `6-tuples.txt`.
     - `Genus 7 generic, part 2.ipynb` find curves with 6 F_2-points. Requires `6-tuples.txt`.
     - `Genus 7 generic, part 3.ipynb` find curves with 7 F_2-points. Requires `6-tuples.txt`.

