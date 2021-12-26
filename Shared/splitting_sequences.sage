# This file accompanies the papers "The relative class number problem for function fields, II"
# by Kiran S. Kedlaya.
#
# This file defines Sage functions manipulate splitting sequences for extensions of function fields.

# Compute the place count of a covering curve from a sequence of splittings.

def place_count_from_splittings(s):
    n = len(s)
    ctp = [0 for _ in range(n)]
    for i in range(n):
        for p in s[i]:
            for j in p:
                if (i+1)*j <= n:
                    ctp[(i+1)*j-1] += 1
    return ctp

# Compute possible splitting sequences given the point counts for $C$ and $C'$ and 
# the degree $d$ of the morphism.
#
# - `alternating`: if $d$ is odd, assume the Galois group is contained in $A_d$; 
#   if $d$ is even, assume the quadratic resolvent is constant.
# - `seed`: if specified, only return sequences that extend one in this list. 
#   Assumes all specified sequences are of the same length.
# - `avoid`: a list of splitting types that must be avoided.
# - `force`: a dictionary of splitting types that must be included in the specified degree.

def splitting_sequences(ct1, ct2, d, alternating=False, seed=[[]], avoid=[], force={}):
    n = len(ct1)
    assert n == len(ct2)
    assert len(set(len(s) for s in seed)) == 1
    ct1p = place_count_from_point_count(ct1, n)
    ct2p = place_count_from_point_count(ct2, n)
    parts = Partitions(d)
    parts0 = [k for k in parts if (len(k)-d)%2 == 0]
    parts1 = [k for k in parts if (len(k)-d)%2 == 1]
    l = seed[:]
    for i in range(len(l[0])+1, n+1):
        l2 = []
        for t in l:
            if not alternating:
                tmp = parts
            elif d%2 == 1 or i%2 == 0:
                tmp = parts0
            else:
                tmp = parts1
            tmp = [p for p in tmp if p not in avoid]
            force_len = 0 if i not in force else len(force[i])
            for u in itertools.combinations_with_replacement(tmp, ct1p[i-1]-force_len):
                u2 = u
                if force_len:
                    u2 = u + force[i]
                t2 = t + [u2]
                if place_count_from_splittings(t2) == ct2p[:i]:
                    l2.append(t2)
        l = l2
    return l
    
# Convert a splitting sequence of a covering to the traces of the quadratic resolvent. 
# This assumes that the quadratic resolvent is purely geometric.

def splittings_to_quadratic_resolvent(ct1, s):
    n = max(len(ct1), len(s))
    d = max(sum(p) for i in range(n) for p in s[i])
    s2 = []
    for t in s:
        s2.append([])
        for p in t:
            if sum(p) == d:
                s2[-1].append([1,1] if (sum(p)-len(p))%2==0 else [2])
            elif d == 3 and p == [1, 1]:
                s2[-1].append([1])
            elif d == 3 and p == [1]:
                s2[-1].append([1,1])
            else:
                raise ValueError("Ramification not supported")
        s2[-1] = tuple(s2[-1])
    ctp = place_count_from_splittings(s2)
    ct = point_count_from_place_count(ctp, n)
    return tuple(ct1[i] - ct[i] for i in range(n))

