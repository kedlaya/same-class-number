# Given a vector space $V$, a subspace $W$ of $V$, and a positive integer $n$, iterate over all the subspaces
# of $V$ containing $W$ with codimension $n$.
# If `basis` is True, return a complementary basis instead of the subspace itself.
def subspaces_containing(V, W, n, basis=False):
    b0 = W.basis()
    b = []
    for v in V.basis():
        if Matrix(b0 + b + [v]).rank() > Matrix(b0 + b).rank():
            b.append(v)
    for W1 in V.subspace(b).subspaces(n):
       if basis:
          yield W1.basis()
       else:
          yield V.subspace(b0 + W1.basis())

# Iterate over solutions of the equation mat*w = v. If sub is specified, it specifies a subspace
# of the kernel, specified as a function of args (so that its execution can be delayed until we
# check whether there are any solutions at all), and we only give representative solutions modulo this subspace.
def solve_right_iterator(mat, v, sub=None, args=None):
    try:
        w = mat.solve_right(v)
        V = mat.right_kernel()
        if sub is None:
            for v1 in V:
                yield v1+w
        else:
            V1 = V.quotient(sub(args))
            for v1 in V1:
                yield V1.lift(v1)+w
    except ValueError:
        pass

