# This file accompanies the papers "The relative class number problem for function fields, I" and
# "The relative class number problem for function fields, II" by Kiran S. Kedlaya.
#
# This file defines Sage functions to construct cyclic covers of function fields 
# with prescribed ramification and to identify cases of relative class number 1 among these.
# These functions depend on Magma via the pexpect interface.

import itertools

# Convert a Magma polynomial to a list of coefficients.
def magma_poly_list(f):
    return [i for i in f.Coefficients()]

# Compute cyclic covers of a Magma function field with fixed ramification divisor.
def cyclic_covers_by_ramification(F, d, M, q=2, delta=0):
    R,mR = M.RayClassGroup(nvals=2)
    gens = [i for i in R.Generators()]
    U1 = R.sub([d*gens[i] for i in range(len(gens))])
    g = Integer(F.Genus())
    g1 = g + (d-1)*(g-1) + delta//(1 if q%2==0 else 2)
    for j in itertools.product(range(d), repeat=len(gens)):
        if all(i==0 for i in j):
            continue
        i = min(i1 for i1 in range(len(gens)) if j[i1] != 0)
        if j[i] > 1 or gcd(d, gcd(j)) > 1:
            continue
        U2 = R.sub([gens[i1] - j[i1]*gens[i] for i1 in range(len(gens))])
        U3 = R.sub([U1, U2])
        A = M.AbelianExtension(U3)
        F1 = magma.FunctionField(A)
        if Integer(F1.Genus()) == g1:
            F1.AssignNames('w')
            yield F1

# Compute cyclic covers of a Magma function field of a fixed degree.
# If delta = 0, the covers must be everywhere unramified.
# If delta > 0, the covers have delta geometric ramification points on the base;
# this assumes delta <= 2.
def cyclic_covers(F, d, delta=0, q=2):
    ans = []
    places1 = [i for i in F.Places(1)]
    num_places1 = len(places1)
    places2 = [i for i in F.Places(2)]
    num_places2 = len(places2)
    z = F.DivisorGroup().Identity()
    if delta == 0:
        m = [z]
    elif delta == 1 and d == 2: #Wild ramification
        m = [z+2*i for i in places1]
    elif delta == 2 and (q%2 == 0 and d == 3) or (q%2 != 0 and d == 2): #Tame ramification
        m = [z+i+j for (i,j) in itertools.combinations_with_replacement(places1, 2)] + \
            [z+i for i in places2]
    elif delta == 2 and d == 2: #Wild ramification
        m = [z+2*i+2*j for (i,j) in itertools.combinations_with_replacement(places1, 2)] + \
            [z+2*i for i in places2]
    elif delta == 1 and d > 3: # This case cannot occur.
        m = []
    else:
        raise ValueError
    for M in m:
        for F1 in cyclic_covers_by_ramification(F, d, M, q=q, delta=delta):
            yield F1

# Exhaust over abelian extensions of a Magma function field of a fixed prime degree to find 
# instances of relative class number 1.

def match_weil_poly(F, d, delta=0, verbose=False):
    ans = []
    for F1 in cyclic_covers(F, d, delta=delta):
        if verbose:
            print(F1.ZetaFunction().Numerator() // F.ZetaFunction().Numerator())
        if F1.ClassNumber() == F.ClassNumber():
            ans.append((F, F1))
    return(ans)

# Reduce a list of Magma function fields to a list of isomorphism class representatives.
def isomorphism_class_reps(l):
    d = defaultdict(list)
    # Hash by point counts.
    for F in l:
        g = Integer(F.Genus())
        t = tuple(Integer(F.NumberOfPlacesOfDegreeOneECF(i)) for i in range(1, g+1))
        d[t].append(F)
    ans = []
    for t in d:
        tmp = []
        for F in d[t]:
            if not any(str(F.IsIsomorphic(F1)) == "true" for F1 in tmp):
                tmp.append(F)
        ans += tmp
    return ans

