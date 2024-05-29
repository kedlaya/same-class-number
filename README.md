# same-class-number

This paper includes code related to the papers "The relative class number one problem for function fields, I, II, III" by Kiran S. Kedlaya. Most of the code is packaged in Jupyter notebooks; these are intended to be run using SageMath (tested using version 9.7). When indicated, there is also an external dependency on Magma (tested using version 2.27-1).

In the subdirectory "Shared":

- `cyclic_covers.sage`: Sage subroutines for finding cyclic covers of function fields. Requires Magma.
- `weil_poly_utils.sage`: Sage subroutines for managing Weil polynomials. This includes converting back and forth between Weil polynomials, point counts, Frobenius traces, and LMFDB labels for isogeny classes of abelian varieties; it also includes calling Sage's Weil polynomials iterator under specific conditions (e.g., with specified Frobenius traces).
- `polys.xlsx`: Excel spreadsheet listing candidate Weil polynomials for purely geometric extensions of function fields over F_2, and indicating which of these correspond to cyclic covers (when known). This is generated and modified by various notebooks.
- `Curves by genus.ipynb`: Process data files from Xarles (see below) and Dragutinović (see below) to build a table of curves of genus up to 5 over F\_2 sorted by their Weil polynomials. Requires Magma. Creates a temporary file `curves.txt`.
- `function_fields.sage`: Loads `curves.txt` and provides an access function to create Magma function fields with a given zeta function.
- `ListCurvesGenus4p2.txt`: data file from [Xarles](https://github.com/XavierXarles/Censusforgenus4curvesoverF2).
- `HyperellipticCurvesData.txt`, `TrigonalCurvesWithAutomorphisms.txt`, `Complete_Intersections.txt`, `Pts_Count_Complete_Intersections.txt`: data files from [Dragutinović](https://github.com/DusanDragutinovic/MT_Curves).

In the subdirectory "Part I":

- `Linear programming.ipynb`: Compute bounds on numbers of rational points on curves over finite fields using the linear programming method.
- `Relative class number 1 for q>2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F_q for q>2, and search for cyclic covers among these. Requires Magma.
- `Weil polynomial bound for q=2.ipynb`: Compute candidate Weil polynomials for purely geometric extensions of function fields over F\_2. Requires `curves.txt`. Updates `polys.xlsx`.
- `Cyclic covers for q=2.ipynb`: Search for cyclic covers among candidate Weil polynomials listed in `polys.xlsx`. Requires Magma and `curves.txt`. Requires and updates `polys.xlsx`.

In the subdirectory "Part II":

- `splitting_sequences.sage`: Sage subroutines for computing splitting sequences for finite extensions of function fields.
- `Non-Galois covers for q=2.ipynb`: Perform various computations needed to rule out the existence of noncyclic covers over F_2. Requires Magma, `curves.txt`, and `polys.xlsx`. 
 
In the subdirectory "Part III":

- `orbits.sage`: Sage subroutines for computing orbit representatives for the action of a finite group on subsets of a finite set. This code is derived from [this repository](https://github.com/kedlaya/orbits).
- `linalg.sage`: Sage subroutines for linear algebra.
- `preamble.sage`: Common declarations for the notebooks in this folder. This includes some initialization code   plus a function `closeout` run at the end to report results.
- A series of Jupyter notebooks covering the Brill-Noether strata in genus 6. These all require Magma and `polys.xlsx`.
  - There is no notebook for hyperelliptic curves because these are handled in the text with no computation at all. Similarly, there is only a short computation needed for bielliptic curves; see below.
  - Trigonal curves are separated by their Maroni invariant: `Genus 6 trigonal Maroni 0.ipynb`, `Genus 6 trigonal Maroni 2.ipynb`.
  - Plane quintics are handled in `Genus 6 plane quintic.ipynb`.
  - Generic curves are handled in a series of three notebooks:
     - `Genus 6 generic, part 1.ipynb` performs the orbit lookup tree calculation. Creates a temporary file `genus6-tuples.txt` used by the other notebooks.
     - `Genus 6 generic, part 2.ipynb` find curves with 5 F\_2-points.
     - `Genus 6 generic, part 3.ipynb` find curves with 4 or 6 F\_2-points.
- A series of Jupyter notebooks covering the Brill-Noether strata in genus 7. These all require Magma and `polys.xlsx`.
  - There is no notebook for hyperelliptic curves because these are handled in the text with no computation at all. Similarly, there is no notebook for trigonal curves of Maroni invariant 3.
  - Trigonal curves of Maroni invariant 1 are handled in `Genus 7 trigonal.ipynb`.
  - Bielliptic curves are handled in `Genus 6 and 7 bielliptic.ipynb`. This notebook also includes a short computation needed to handle bielliptic curves in genus 6.
  - Curves with a g\^2\_6 are separated based on whether the g\^2\_6 is self-adjoint, rational over F\_2, or irrational over F\_2: `Genus 7 self-adjoint g26.ipynb`, `Genus 6 rational g26.ipynb`, `Genus 6 irrational g26.ipynb`.
  - Tetragonal curves (with no g\^2\_6) are handled in `Genus 7 tetragonal, no g26.ipynb`.
  - Generic curves are handled in a series of three notebooks:
     - `Genus 7 generic, part 1.ipynb` performs the orbit lookup tree calculation. Creates a temporary file `6-tuples.txt` used by the other notebooks.
     - `Genus 7 generic, part 2.ipynb` find curves with 6 F\_2-points.
     - `Genus 7 generic, part 3.ipynb` find curves with 7 F\_2-points.

The author hereby consigns all original code and data included in this repository to the public domain.
