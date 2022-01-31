# This file accompanies the papers "The relative class number problem for function fields, I, II, III" by Kiran S. Kedlaya.
#
# This file loads the dictionary `curves_by_poly` specifying all curves of genus up to 5 (and selected curves of 
# genus 6 and 7) from the file "curves.txt". It also defines a Sage function `function_fields_by_weil_poly` to convert
# these curves into Magma function fields of curves with a specified zeta function.

with open('../Shared/curves.txt', 'r') as f:
    s = f.read()
P4.<x0,x1,x2,x3,x4> = GF(2)[]
curves_by_poly = sage_eval(s, locals={'x': x0, 'y': x1, 'x0': x0, 'x1': x1, 'x2': x2, 'x3': x3, 'x4': x4})
P0.<x> = GF(2)[]
for g in range(1, 8):
    for s in curves_by_poly[g]:
        m = curves_by_poly[g][s] if g <= 2 else curves_by_poly[g][s][0]
        for i in range(len(m)):
            (h, f) = m[i]
            if h == 1:
                h = P4(1)
            if f == 1:
                f = P4(1)
            m[i] = (h(x,0,0,0,0), f(x,0,0,0,0))
P3.<x0,x1,x2> = GF(2)[]
for s in curves_by_poly[3]:
    curves_by_poly[3][s][1] = [f(x0,x1,x2,0,0) for f in curves_by_poly[3][s][1]]
for s in curves_by_poly[5]:
    curves_by_poly[5][s][1] = [f(x0,x1,x2,0,0) for f in curves_by_poly[5][s][1]]
P4.<x0,x1,x2,x3> = GF(2)[]
for s in curves_by_poly[4]:
    curves_by_poly[4][s][1] = [tuple(f(x0,x1,x2,x3,0) for f in t) for t in curves_by_poly[4][s][1]]
Pxy.<x,y> = GF(2)[]
for g in range(6, 8):
    for s in curves_by_poly[g]:
        curves_by_poly[g][s][1] = [f(x,y,0,0,0) for f in curves_by_poly[g][s][1]]

def function_fields_by_weil_poly(s):
    F0 = magma.FunctionField(GF(2))
    F0.AssignNames('x')
    S = magma.PolynomialRing(F0)
    S.AssignNames('y')
    l = s.split(".")
    g = Integer(l[0])
    if s in curves_by_poly[g]:
        ans = []
        # Hyperelliptics.
        m = curves_by_poly[g][s] if g <= 2 else curves_by_poly[g][s][0]
        for (u,v) in m:
            yield magma.HyperellipticCurve(v,u).FunctionField().AlgorithmicFunctionField()
        # Plane quartics for g=3.
        if g == 3:
            for u in curves_by_poly[g][s][1]:
                proj = magma.ProjectiveSpace(u.parent())
                yield proj.Scheme([u]).FunctionField().AlgorithmicFunctionField()
        # Genus 4 curves.
        if g == 4:
            for (u,v) in curves_by_poly[g][s][1]:
                proj = magma.ProjectiveSpace(u.parent())
                yield proj.Scheme([u,v]).FunctionField().AlgorithmicFunctionField()
        # Genus 5 curves.
        if g == 5:
            for u in curves_by_poly[g][s][1]:
                proj = magma.ProjectiveSpace(u.parent())
                yield proj.Scheme([u]).FunctionField().AlgorithmicFunctionField()        
            for (u,v,w) in curves_by_poly[g][s][2]:
                proj = magma.ProjectiveSpace(u.parent())        
                yield proj.Scheme([u,v,w]).FunctionField().AlgorithmicFunctionField()
        # Genus 6 and 7 curves, represented as singular plane curves.
        if g == 6 or g == 7:
            for u in curves_by_poly[g][s][1]:
                proj = magma.AffineSpace(u.parent())
                yield proj.Scheme(u).FunctionField().AlgorithmicFunctionField()
