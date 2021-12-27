# Given a group G, construct a random generating sequence. This is likely to be shorter than a generating sequence
# found in some other way, and hence more efficient for constructing Cayley graphs.

def random_generating_sequence(G):
    l = []
    while G.subgroup(l).order() != G.order():
        l.append(G.random_element())
    return l
    
# Given a subspace of $F^n$ represented by an $m \times n$ matrix, compute the subgroup of $GL(n,F)$ stabilizing 
# this subspace for the *right* action. Make sure to transpose if you want a left action!

def vec_stab(M):
    F = M.base_ring()
    m,n = M.dimensions()
    # Conjugate to a standard matrix.
    l = M.rows()
    for i in range(n):
        if i not in M.pivots():
            l.append(vector(F, (1 if j==i else 0 for j in range(n))))
    M1 = Matrix(F, l)
    # Then construct a block matrix group.
    l0 = [block_matrix(2,2,[g.matrix(),0,0,identity_matrix(n-m)], subdivide=False) for g in GL(m, F).gens()] + \
        [block_matrix(2,2,[identity_matrix(m),0,0,g.matrix()], subdivide=False) for g in GL(n-m, F).gens()]
    l0.append(identity_matrix(n))
    l0[-1][m,0] = 1
    G = GL(n, F).subgroup([~M1*g*M1 for g in l0])
    assert all((M*g.matrix()).echelon_form() == M.echelon_form() for g in G.gens())
    return G

# Given a generalized Cayley graph for a group G acting on a set of vertices, return a list of connected
# component representatives not containing any vertices listed in `exclude` or any vertices 
# for which `forbid` returns True, and a dictionary computing these representatives to 
# arbitrary vertices in their components.

def group_retract(G, vertices, edges, exclude=[], forbid=None):
    vertices = list(vertices) + exclude
    Gamma = Graph([vertices, edges], loops=True, format='vertices_and_edges')
    assert all(x in vertices for x in Gamma)
    conn = Gamma.connected_components()
    reps = [l[0] for l in conn]
    if exclude and forbid:
        for M in exclude + reps:
            if (M in exclude) or (M in reps and forbid(M)):
                Gamma.add_edge(exclude[0], M)
    forbidden_verts = []
    for v in exclude:
        forbidden_verts += Gamma.connected_component_containing_vertex(v, sort=False)
    Gamma.delete_vertices(forbidden_verts)
    conn = Gamma.connected_components()
    reps = [l[0] for l in conn]
    d = {M: (M, G(1)) for M in reps}
    queue = reps[:]
    while queue:
        M = queue.pop(0)
        assert M in d
        nb = [M1 for M1 in Gamma.neighbors(M) if M1 not in d]
        for M1 in nb:
            assert M1 not in queue
            queue.append(M1)
            g = Gamma.edge_label(M, M1)
            if (M1, M, g) in edges:
               g = ~g
            elif (M, M1, g) not in edges:
                raise AssertionError
            d[M1] = (d[M][0], g*d[M][1])
    assert all(M in d or M in forbidden_verts for M in vertices)
    return reps, d, forbidden_verts

# The data structure for orbit lookup trees of depth $n$:
# - The tree is a dictionary `tree` indexed by levels $0, \dots, n$.
# - Each level `tree[k]` is a dictionary keyed by a $k$-tuple, identified with nodes of the tree.
# - Each value `tree[k][U]` is a dictionary with the following keys:
#  - `gpel` (for $U$ eligible): a pair $(U', g)$ such that $U'$ is a green node and $g(U') = U$.
#  - `stab` (for $U$ green) the stabilizer of $U$.
#  - `retract` (for $U$ green and $k<n$): a dictionary whose value at $y \in S \setminus U$ (resp $y \in S/U$) is an element $g \in G_U$ such that $U \cup \{g^{-1}(y)\}$ (resp. $\pi_U^{-1}(g^{-1}(y))$) is a red or green node.

# Use an enhanced $n$-orbit tree to identify an orbit representative for the action of the group $G$ on $k$-tuples of Lagrangian subspaces. 
#
def orbit_rep_from_tree(G, tree, mats, apply_group_elem, find_green=True):
    n = len(mats)
    if n not in tree:
        raise ValueError("Tree not computed")
    if n == 0:
        return mats, G(1)
    mats0 = mats[:-1]
    if mats0 in tree[n-1] and 'gpel' not in tree[n-1][mats0]: # Truncation is an ineligible node
        return None, None
    if mats0 in tree[n-1] and 'stab' in tree[n-1][mats0]: # Truncation is a green node
        assert 'retract' in tree[n-1][mats0]
        g0 = G(1)
        y = mats[-1]
        assert y not in mats0
    else: # Truncation needs to be resolved
        mats0, g0 = orbit_rep_from_tree(G, tree, mats[:-1], apply_group_elem)
        if mats0 is None:
            return None, None
        assert 'gpel' in tree[n-1][mats0]
        assert 'retract' in tree[n-1][mats0]
        y = apply_group_elem(~g0, mats[-1])
    z, g1 = tree[n-1][mats0]['retract'][y]
    assert z not in mats0
    mats1 = mats0 + (z,)
    if not find_green:
        return mats1, g0*g1
    if 'gpel' not in tree[n][mats1]:
        return None, None
    mats2, g2 = tree[n][mats1]['gpel']
    assert 'gpel' in tree[n][mats2]
    return mats2, g0*g1*g2

# Given an orbit lookup tree at depth $n$ (for the action of a finite group $G$ on a finite set $S$), extend it in place
# to depth $n+1$. For $n=0$, pass for `tree` an empty dict and it will be initialized correctly.
#
# The argument `methods` is a dictionary with the following keys:
# - `apply_group_elem`: given a pair $(g, x) \in G \times S$, returns $g(x)$.
# - `stabilizer`: given $x \in S$, returns a group whose intersection with $G$ (in some ambient group) is $G_x$.
# - `forbid` (optional): given a tuple $(x_1,\dots,x_k)$, return True if the underlying subset $\{x_1,\dots,x_k\}$ is forbidden. It is assumed that this function is symmetric in the input tuple. If some of these checks are time-consuming them, only run them when the optional argument `easy` is True.

def extend_orbit_tree(G, S, tree, methods, verbose=True, terminate=False):
    apply_group_elem = methods['apply_group_elem']
    stabilizer = methods['stabilizer']
    forbid = methods['forbid'] if 'forbid' in methods else (lambda x, easy=False: False)
    if not tree: # Initialize
        S0 = tuple()
        tree[0] = {S0: {'gpel': (S0, G(1)), 'stab': []}}
    n = max(tree.keys())
    if verbose:
        print("Current level: {}".format(n))
    tree[n+1] = {}
    for mats in tree[n]:
        if 'stab' in tree[n][mats]: # Green node
            # Compute the stabilizer of mats (deferred from the previous level).
            if n == 0:
                G0 = G
            else:
                parent = mats[:-1]
                endgen = mats[-1]
                G0 = tree[n-1][parent]['stab'].intersection(stabilizer(endgen))
            G1 = G.subgroup(list(G0.gens()) + tree[n][mats]['stab'])
            vertices = [M for M in S if M not in mats]
            if G1.order() == 1: # Early abort
                if verbose:
                    print("Trivial stabilizer")
                tmp = vertices
                d = {M: (M, G(1)) for M in vertices}
            else: # Construct the group retract under this green node.
                gens = random_generating_sequence(G1)
                assert all(apply_group_elem(g, M) in mats for g in gens for M in mats)
                if verbose:
                    print("Number of generators: {}".format(len(gens)))
                G1 = G.subgroup(gens)
                edges = [(M, apply_group_elem(g, M), g) for M in vertices for g in gens]
                tmp, d, _ = group_retract(G, vertices, edges)
            tree[n][mats]['stab'] = G1
            tree[n][mats]['retract'] = d
            for M in tmp:
                if M in mats:
                    raise ValueError("Found repeated entry in tuple")
                mats1 = mats + (M,)
                tree[n+1][mats1] = {}
    if verbose:
        print("Number of new nodes: {}".format(len(tree[n+1])))
    dummy_vertex = tuple()
    edges = []
    exclude = [dummy_vertex]
    for mats in tree[n+1]:
        if forbid(mats, easy=True):
            exclude.append(mats)
        else:
            tmp = [tuple(mats[n if i==j else j if i==n else i] for i in range(n+1)) for j in range(n)]
            for i in tmp:
                mats1, g1 = orbit_rep_from_tree(G, tree, i, apply_group_elem, find_green=False)
                if mats1 is None:
                    exclude.append(mats)
                else:
                    edges.append((mats1, mats, g1))
    if verbose:
        print("Edges computed")
    tmp, d, forbidden_verts = group_retract(G, tree[n+1].keys(), edges, exclude, forbid)
    if verbose:
        print("Retract computed")
    for mats in d:
        tree[n+1][mats]['gpel'] = d[mats]
    # Mark green nodes.
    for mats in tmp:
        tree[n+1][mats]['stab'] = []
    assert all('stab' in tree[n+1][tree[n+1][mats]['gpel'][0]] 
                   for mats in tree[n+1] if 'gpel' in tree[n+1][mats])
    # Defer the computation of stabilizers, recording some elements read off of the graph.
    if terminate:
        if verbose:
            print("Stabilizer generators not computed")
    else:
        for e in edges:
            if e[0] not in forbidden_verts:
                mats1, g1 = d[e[0]]
                mats2, g2 = d[e[1]]
                assert mats1 == mats2
                g = ~g2*e[2]*g1
                if g != G(1):
                    tree[n+1][mats1]['stab'].append(g)
        if verbose:
            print("Stabilizer generators found")
    if verbose:
        print("Number of new green nodes: {}".format(sum(1 for mats in tree[n+1] 
                                                     if 'stab' in tree[n+1][mats])))
        print("New level: {}".format(max(tree.keys())))

# Build an orbit lookup tree to depth n. By default, we do not record stabilizer generators at depth n,
# so the result cannot be used to extend to depth n+1.

def build_orbit_tree(G, S, n, methods, verbose=True, terminate=True):    
    # Verify that each generator of G defines a permutation of S.
    apply_group_elem = methods['apply_group_elem']
    for g in G.gens():
        Sg = [apply_group_elem(g, x) for x in S]
        assert all(x in S for x in Sg) and all(x in Sg for x in S)
    tree = {}
    for i in range(n):
        extend_orbit_tree(G, S, tree, methods, verbose=verbose, terminate=(terminate and (i == n-1)))
    return tree

# Return a list of green nodes at depth $n$.
def green_nodes(tree, n):
    return [mats for mats in tree[n] if 'stab' in tree[n][mats]]
