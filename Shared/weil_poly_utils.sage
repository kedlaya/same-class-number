# This file accompanies the papers "The relative class number problem for function fields, I" and
# "The relative class number problem for function fields, II" by Kiran S. Kedlaya.
#
# This file defines Sage functions to manipulate Weil polynomials, their point counts, and their
# Frobenius traces.

import itertools
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials

# Given the Weil polynomial of an abelian variety, return the first n Frobenius traces.
def trace_from_weil_poly(u, n):
    Q.<t> = PowerSeriesRing(QQ)
    v = u.reverse()(t).log()
    l = v.list()
    while len(l) <= n+1:
        l.append(0)
    return [-l[i]*i for i in range(1, n+1)]
 
# Given the Weil polynomial of a curve over F_q, return the first n point counts.
def point_count_from_weil_poly(u, n, q=2):
    tmp = trace_from_weil_poly(u, n)
    return [q^i+1-tmp[i-1] for i in range(1, n+1)]

# Given the first n point counts of a curve over F_q, return the first n place counts.
def place_count_from_point_count(s, n):
    return [sum(moebius(j)*s[i//j-1] for j in divisors(i))//i for i in range(1, n+1)]
    
# Given the Weil polynomial of a curve over F_q, return the first n place counts.
def place_count_from_weil_poly(u, n, q=2):
    return place_count_from_point_count(point_count_from_weil_poly(u, n, q=q), n)
    
# Given the Weil polynomial of a curve over F_q, return True if its first n place counts are nonnegative.
def check_curve_positivity(u, n, q=2):
    s = place_count_from_weil_poly(u, n, q=q)
    return all(_ >= 0 for _ in s)

# Given the first n place counts of a curve over F_q, return the first n point counts.
def point_count_from_place_count(s, n):
    return [sum(s[j-1]*j for j in divisors(i)) for i in range(1,n+1)]
    
# Given point counts of a curve of a given genus, return its Weil polynomial.
def weil_poly_from_point_count(l, d, q=2):
    P.<T> = QQ[]
    Q.<t> = PowerSeriesRing(QQ)
    u = sum(l[i-1]*t^i/i for i in range(1,d+1))
    v = exp(u)*(1-t)*(1-q*t)
    l2 = [v[i] for i in range(d+1)]
    for i in range(1, d+1):
        l2.append(2^i*v[d-i])
    return P(l2).reverse()

# Identify all Weil polynomials with specified initial Frobenius traces.
def weil_polys_from_traces(P, d, q, l):
    Q.<t> = PowerSeriesRing(QQ)
    u = sum(-l[i-1]*t^i/i for i in range(1,len(l)+1))
    v = exp(u)
    l2 = v.polynomial().list()[:len(l)+1]
    l2 = [ZZ(i) for i in l2]
    if len(l2) > d//2 + 1:
        l2 = l2[:d//2 + 1]
    for w in WeilPolynomials(d=d, q=q, lead=l2):
        if w.is_weil_polynomial() and trace_from_weil_poly(w, len(l)) == list(l):
            yield P(w)
#    l3 = P.weil_polynomials(d=d, q=q, lead=l2)
#    return [w for w in l3 if w.is_weil_polynomial() and trace_from_weil_poly(w, len(l)) == list(l)]

# Identify all Weil polynomials with specified initial point counts.
def weil_polys_from_point_counts(P, d, q, l):
    return weil_polys_from_traces(P, d, q, [q^i+1-l[i-1] for i in range(1, len(l)+1)])

# Return the LMFDB label corresponding to a given Weil polynomial.
def label_from_weil_poly(f):
    q = f.trace_polynomial()[2]
    d = f.degree()
    l = list(f)
    g = d // 2
    s = '{}.{}.'.format(g, q)
    for i in range(g):
        if l[-i-2] <= 0:
            s += 'a'
        t = Integer(abs(l[-i-2])).digits(26)
        for j in reversed(t):
            s += chr(97+j)
        if i<g-1:
            s += '_'
    return s
    
# Return the Weil polynomial corresponding to an LMFDB label.
def weil_poly_from_label(P, s):
    l = s.split(".")
    g = int(l[0])
    q = int(l[1])
    l2 = l[2].split("_")
    l3 = [1]
    for s1 in l2:
        sg = -1 if s1[0] == 'a' else 1
        t = [ord(j)-97 for j in s1]
        t.reverse()
        l3.append(sg*ZZ(t, base=26))
    l4 = l3 + [q^i*l3[-i-1] for i in range(1,len(l3))]
    l4.reverse()
    return P(l4)
    
# Compute the "modified reduced resultant" of two real Weil polynomials in the sense of Howe-Lauter.
# This function is taken from LMFDB and is originally due to Everett Howe (ported from Magma).
def modified_reduced_resultant(h1, h2, q=2):
    """
    Suppose ``h1`` and ``h2`` are real Weil `q`-polynomials (assumed to be coprime to one
    another), with associated Frobenius elements ``pi1`` and ``pi2`` in the centers of
    their associated endomorphism rings ``R1`` and ``R2``.  Let ``pi`` be the element
    ``(pi1, pi2)`` of ``R1 x R2``, and let ``pibar = q/pi``.  The *modified reduced
    resultant* of ``h1`` and ``h2`` is the smallest positive integer `n` such that ``(0, n)``
    lies in the subring ``Z[pi,pibar]`` of ``R1 x R2``.
    Usually, this `n` is simply the reduced resultant of ``h1`` and ``h2``, but if the
    product ``h1*h2`` is divisible by ``x^2 - 4*q``, we can sometimes divide the reduced
    resultant by 2.
    """
    T = h1.parent().gen()
    d, a1, a2 = xgcd(h1.change_ring(QQ), h2.change_ring(QQ))
    if d.degree() > 0:
        return 0
    n = lcm(c.denominator() for c in a1)
    h = h1*h2
    H, rem = h.quo_rem(T^2-4*q)
    if rem == 0:
        splitelt = n*a1*h1 # 0 mod h1, n mod h2
        otherelt = splitelt + T*H
        g = gcd([ZZ(c) for c in otherelt])
        if g % 2 == 0:
            return n // 2
    return n
    
# Given a Weil polynomial, return True if the Serre-Howe-Lauter resultant criterion indicates that there is
# no Jacobian over F_q with this Weil polynomial.
# The resultant 1 case is based on code from LMFDB due to Everett Howe (ported from Magma).
# We also implement a form of the resultant 2 case, using the Castelnuovo-Severi inequality and
# the Deuring-Shafarevich formula.
# If blocklist is specified, it should be a dictionary in which blocklist[d] is a list of Weil polynomials
# of degree 2*d which are known not to occur for Jacobians.
def _nojac_serre(pol, q=2, alert=False, blocklist=None):
    pol0 = pol.trace_polynomial()[0]
    irred_factors = [f for (f, _) in pol0.factor()]
    n = len(irred_factors)
    res2list = []
    for i in range(2^(n-1)-1):
        h0 = h1 = 1
        for b, factor in zip((2^(n-1)+i).bits(), irred_factors):
            if b == 0:
                h0 *= factor
            else:
                h1 *= factor
        res = modified_reduced_resultant(h0, h1, q=q)
        if res == 1:
            return True
        if res == 2:
            h0full = prod(f^e for (f,e) in pol0.factor() if h0%f == 0)
            h1full = pol0 // h0full
            for h in (h0full, h1full):
                delta = pol0.degree() - 2*h.degree() + 1
                rank0 = max(pol0.degree()-i for i in range(pol0.degree()+1) if pol0[i]%2)
                rank1 = max(h.degree()-i for i in range(h.degree()+1) if h[i]%2)
                h2 = h.reciprocal_transform(q=2)
                if delta >= 0 and check_curve_positivity(h2, h.degree()) and \
                    (q%2 or (rank0 - 2*rank1 + 1 in range(delta+1))) and \
                    (not blocklist or h2 not in blocklist[h.degree()]):
                    # Check the Castelnuovo-Severi inequality.
                    res2list.append(h.degree())
                    if len(res2list) >= 2:
                        res2list.sort()
                        if 2*res2list[0] + 2*res2list[1] < pol.degree():
                            if alert:
                                print("Alert: Castelnuovo-Severi inequality triggered")
                            return True
                    if alert:
                        print("Alert: {} {}".format(pol, h2))
                    break
            else:
                return True
        if alert:
            print(res, h0, h1)
        if alert and res <= 10:
            print("Alert: {} {} {} {}".format(res, pol, h0.reciprocal_transform(q=2), \
                                                        h1.reciprocal_transform(q=2)))
   

