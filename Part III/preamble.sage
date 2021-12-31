## Common preamble for genus 6 and 7 enumerations

# Import packages

import itertools, pandas
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
        
# Extract the relevant point counts in genus 6 and genus 7.
# The resulting tuples are in the lists `targets6` and `targets7`.

targets6 = []
targets7 = []
for i in range(len(df)):
    r = df.iloc[i]
    d = r["d"]
    g = r["g"]
    g1 = r["g'"]
    if (d,g,g1) != (2,6,11) and (d,g,g1) != (2,7,13):
        continue
    counts1 = tuple(eval(r["Counts of C"]))
    if g == 6:
        targets6.append(counts1[:6])
    elif g == 7:
        targets7.append(counts1[:7])
assert all(t[0] in [4,5,6] for t in targets6)
assert all(t[0] in [6,7] for t in targets7)


