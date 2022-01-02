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

# Record the current time.

timestamp = time.time()

# Declare a function to close out a given notebook.
# The input is assumed to be a list of Magma curves if fields is False (default), or function fields if fields is True.

def report_time():
    seconds = ceil(time.time() - timestamp) 
    if seconds < 300:
        print("Total time: {} seconds".format(seconds))
    else:
        minutes = ceil(seconds//60)
        print("Total time: {} minutes".format(minutes))

def closeout(l = [], fields=False):
    if not l:
        print("No curves found in this case!")
        report_time()
        return
    # Convert to function fields.
    if fields:
        l2 = [F.AlgorithmicFunctionField() for F in l]
    else:
        l2 = [C.FunctionField().AlgorithmicFunctionField() for C in l]
    # Infer the genus.
    genus = Integer(l2[0].Genus())
    assert all(F.Genus() == genus for F in l2)
    # Identify distinct isomorphism classes.
    l3 = isomorphism_class_reps(l2)
    print("Number of curves found: {}".format(len(l3)))
    # Search for extensions with relative class number 1.
    l4 = []
    for F1 in l3:
        l4 += match_weil_poly(F1, 2)
    if not l4:
        print("No covers found in this case!")
        report_time()
        return
    print("Number of covers found: {}".format(len(l4)))
    # Display our results.
    for (F1, F2) in l4:
        F1.AssignNames("y")
        t = tuple(F1.NumberOfPlacesOfDegreeOneECF(i) for i in range(1, genus+1))
        print(t, F1.RationalExtensionRepresentation())
    # Convert into numerical data in preparation for updating the spreadsheet.
    l5 = []
    Q.<T> = QQ[]
    for (F1, F2) in l4:
        u1 = Q(magma_poly_list(F1.ZetaFunction().Numerator())).reverse()
        v1 = Q(magma_poly_list(F2.ZetaFunction().Numerator())).reverse()
        l5.append((2, genus, 2*genus-1, point_count_from_weil_poly(u1, 13), point_count_from_weil_poly(v1, 13)))
    # Write the results back to the spreadsheet. In the process, we check that every covering we found gives a 
    # pair of Weil polynomials from our original list.
    used_pols = []
    for i in range(len(df)):
        r = df.iloc[i]
        d = r["d"]
        g = r["g"]
        g1 = r["g'"]
        ct1 = eval(r["Counts of C"])
        ct2 = eval(r["Counts of C'"])
        tmp = (d, g, g1, ct1, ct2)
        if tmp in l5:
            df.loc[i, "Cyclic"] = "Yes"
            used_pols.append(tmp)
        if df.loc[i, "Cyclic"] == "Unknown":
            df.loc[i, "Cyclic"] = "No"
    assert all(t in used_pols for t in l5)
    df.to_excel('../Shared/polys.xlsx', sheet_name='Weil polynomials', merge_cells=True, freeze_panes=(int(1),int(1)))
    # Announce that we're done, and report the timing.
    print("All covers recorded!")
    report_time(timestamp)

