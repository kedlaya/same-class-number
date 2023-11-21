## Common preamble for genus 5, 6, and 7 enumerations

# Import packages

import itertools, pandas, time
from collections import defaultdict

# Load other Sage files

load("orbits.sage") # Group orbits
load("linalg.sage") # Auxiliary linear algebra
load("../Shared/cyclic_covers.sage") # Cyclic covers of function fields
load("../Shared/weil_poly_utils.sage") # Utility functions for Weil polynomials

# Retrieve the table of point counts from the Excel spreadsheets.

df = pandas.read_excel('../Shared/polys.xlsx')
for s in df:
    if 'Unnamed' in s:
        del df[s]
        
# Extract the relevant point counts in genus 5, 6, and 7.
# The resulting tuples are in the lists `targets5`, `targets6`, and `targets7`.

targets5 = []
targets6 = []
targets7 = []
for i in range(len(df)):
    r = df.iloc[i]
    d = r["d"]
    g = r["g"]
    g1 = r["g'"]
    if (d,g,g1) not in [(2,5,9), (2,6,11), (2,7,13)]:
        continue
    counts1 = tuple(eval(r["Counts of C"]))
    if g == 5:
        targets5.append(counts1[:5])
    elif g == 6:
        targets6.append(counts1[:6])
    elif g == 7:
        targets7.append(counts1[:7])
assert all(t[0] in [2,3,4,5,6] for t in targets5)
assert all(t[0] in [4,5,6] for t in targets6)
assert all(t[0] in [6,7] for t in targets7)

targets_dict = {5: targets5, 6: targets6, 7: targets7}

# Record the current time.

timestamp = time.time()

# Declare a function to close out a given notebook.
# The input is assumed to be a dictionary `curves` indexed by tuples of positive integers.
# The value `curves[s]` is assumed to be a list of tuples of generators of subschemes of the Magma projective scheme
# `X` (if `X` is specified) or Magma function fields (if `X` is omitted).
# If `genus` is specified, verify that the curves all have the specified genus.

def report_time():
    seconds = ceil(time.time() - timestamp) 
    if seconds < 300:
        print("Total time: {} seconds".format(seconds))
    else:
        minutes = ceil(seconds//60)
        print("Total time: {} minutes ({} seconds)".format(minutes,seconds))

def closeout(curves=None, X=None, genus=None):
    if not curves:
        print("No curves found in this case!")
        report_time()
        return
    # Pick out cases with matching point counts.
    l = []
    Q.<T> = QQ[]
    for s in curves:
        for gens in curves[s]:
            if not X: # In this case, assume gens is already a function field.
                F = gens
            else:
                Y = X.Scheme(gens)
                # Only keep Y if it is integral of dimension 1.
                if Y.Dimension() > 1 or str(Y.IsIrreducible()) == "false" or str(Y.IsReduced()) == "false":
                    continue
                C = Y.Curve()
                F0 = C.FunctionField()
                # Convert F0 into Magma's preferred internal representation.
                F = F0.AlgorithmicFunctionField()
            g = Integer(F.Genus())
            if genus and (genus != g):
                continue
            ct = tuple(Integer(F.NumberOfPlacesOfDegreeOneECF(i)) for i in range(1, g+1))
            # Cross-check the point counts we already computed.
            assert ct[:len(s)] == s
            if ct in targets_dict[g]:
                l.append(F)
    print("Number of curves found: {}".format(len(l)))
    # Identify distinct isomorphism classes.
    l2 = isomorphism_class_reps(l)
    print("Number of isomorphism classes found: {}".format(len(l2)))
    # Search for extensions with relative class number 1.
    l3 = []
    for F1 in l2:
        l3 += match_weil_poly(F1, 2)
    if not l3:
        print("No covers found in this case!")
        report_time()
        return
    print("Number of covers found: {}".format(len(l3)))
    # Display our results.
    for (F1, F2) in l3:
        F1.AssignNames("y")
        g = Integer(F1.Genus())
        t = tuple(F1.NumberOfPlacesOfDegreeOneECF(i) for i in range(1, g+1))
        print(t, F1.RationalExtensionRepresentation().DefiningPolynomial(), 
            F1.AutomorphismGroup().IdentifyGroup())
    # Convert into numerical data in preparation for updating the spreadsheet.
    l4 = []
    for (F1, F2) in l3:
        g = Integer(F1.Genus())
        u1 = Q(magma_poly_list(F1.ZetaFunction().Numerator())).reverse()
        v1 = Q(magma_poly_list(F2.ZetaFunction().Numerator())).reverse()
        l4.append((2, g, 2*g-1, point_count_from_weil_poly(u1, 13), point_count_from_weil_poly(v1, 13)))
    # Verify that every cover we found is already accounted for in the spreadsheet.
    for i in range(len(df)):
        r = df.iloc[i]
        d = r["d"]
        g = r["g"]
        g1 = r["g'"]
        ct1 = eval(r["Counts of C"])
        ct2 = eval(r["Counts of C'"])
        tmp = (d, g, g1, ct1, ct2)
        if tmp in l4:
            assert df.loc[i, "Cyclic"] == "Yes"
    # Announce that we're done, and report the timing.
    print("All covers recorded!")
    report_time()
