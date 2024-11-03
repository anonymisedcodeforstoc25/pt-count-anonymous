from itertools import chain, combinations
from typing import Dict, Iterable, List
from sage.all import Permutations, Permutation, Matrix, QQ, Composition, \
    Compositions, Graph, compose
from sage.graphs.trees import TreeIterator


PermInt = Dict[Permutation, int]


def main():
    ''' Compute pattern-trees matrices mentioned in the paper. '''

    print('****** O(n^2 polylog n) for k <= 7 ******')
    print('Number of permutations: 5913 = 1! + 2! + 3! + 4! + 5! + 6! + 7!')
    # Save time by restricting tree order to 4, which happens to be sufficient.
    show_rank(max_k=7, max_v_size=2, max_num_verts=4, special5=False)

    print('\n****** O(n^(7/4) polylog n) for k <= 5 ******')
    print('Number of permutations: 153 = 1! + 2! + 3! + 4! + 5!')
    show_rank(max_k=5, max_v_size=1, max_num_verts=5, special5=True)


def show_rank(max_k: int, max_v_size: int, max_num_verts: int, special5: bool):
    ''' Compute and print the rank corresponding to the given bounds. '''

    # Get all the needed rows (with the special k<=5 ones if needed)
    # and feed them to the next function.
    # Each row is represented sparsely as a PermInt: a mapping of permutations
    # to integers (that is, a linear combination of patterns).
    rows: Iterable[PermInt] = gen_rows(max_k, max_num_verts, max_v_size)
    if special5:
        rows = chain(rows, gen_special5_rows())
    show_rank_from_rows(max_k, rows)


def show_rank_from_rows(max_k: int, rows: Iterable[PermInt]):
    ''' Compute and print the rank of the matrix containing the given rows. '''

    # The matrix has |S_{<=max_k}| columns, indexed by perms_idx
    # (a mapping from permutations of size <= max_k to integers).
    # The number of rows grows as needed.
    perms_idx = enum_perms(max_k)
    num_cols = len(perms_idx)
    mat = Matrix(QQ, 1, num_cols, sparse=True)

    # Generate sparse rows and add them to the matrix.
    for row_idx, row in enumerate(rows):
        # Check whether more space needs to be allocated.
        if row_idx >= mat.nrows():
            print(' ' * 50, '\rGenerated %s rows so far.' % mat.nrows(),
                  'Rank: ', end='', flush=True)
            cur_rank = mat.rank()
            print(cur_rank, end='\r', flush=True)

            # If the matrix has full rank, we can exit early.
            if cur_rank == num_cols:
                break

            # Double the amount of available rows by padding with zeros.
            mat = mat.stack(Matrix.zero(QQ, row_idx, num_cols, sparse=True))

        # Add the current row. It is a mapping of permutations to
        # integers, representing a linear combination of pattern counts.
        for p, val in row.items():
            mat[row_idx, perms_idx[p]] += val

    print(' ' * 50, end='\r')
    print('Total %s rows. Rank:' % mat.nrows(), end=' ', flush=True)
    print(mat.rank())


def gen_rows(max_k: int, max_num_verts: int,
             max_v_size: int) -> Iterable[PermInt]:
    ''' Generate rows given bounds on the pattern-tree parameters. '''

    # A direct implementation would enumerate over every pattern-tree
    # and then compute its linear combination from linear extensions,
    # as explained in Section 3. This only needs to be computed once
    # to obtain the full-rank matrix.
    #
    # For speed, we do this in a different order. We condition
    # on the size of the pattern ("k"), the order of the tree ("num_verts"),
    # the unlabeled graph structure ("tree"), and the number of pattern
    # points in each vertex ("comp" - a composition is an ordered
    # list of integers that sum to k).
    #
    # This is faster because once we condition on these parameters,
    # the relative order of x- and y- coordinates of the pattern points
    # fully determine the pattern-tree's labels - we don't have to search
    # in order to match patterns to pattern-trees.
    #
    # Those labels are computed in get_pt_id. Note that gen_ptc_rows
    # is called separately for different values of k and tree structures,
    # so the IDs generated are consistent inside each call.
    return chain.from_iterable(
        gen_ptc_rows(k, tree, comp)
        for k in range(1, max_k + 1)
        for num_verts in range(1, max_num_verts + 1)
        for tree in TreeIterator(num_verts)
        for comp in Compositions(k, length=num_verts, max_part=max_v_size)
    )


def gen_ptc_rows(k: int, tree: Graph, comp: Composition) -> Iterable[PermInt]:
    ''' Generate rows for a given tree structure and vertex sizes
        (a composition of k) '''

    # A mapping of pattern-tree IDs (explained below) to linear combinations
    # of pattern counts.
    res: Dict[int, PermInt] = {}

    # The pattern points are indexed 0,1,...,k-1. v_idxs is a map
    # from every pattern-tree vertex to the pattern points it contains.
    v_idxs = get_v_idxs(tree, comp)

    # Go over every relative ordering of the k points in x and y.
    for perm_x in Permutations(k):
        # We do not count relative x-orderings in which points
        # inside a vertex are not sorted (see formalism in Section 3).
        if not allow_perm_x(v_idxs, perm_x):
            continue
        for perm_y in Permutations(k):
            # The pattern formed by this ordering is the composition
            # of y-ordering on the inverse of x-ordering (see Lemma 3.5).
            p = perm_y.left_action_product(perm_x.inverse())

            # Conditioned on the given tree structure and vertex sizes,
            # these orderings determine all labels, so they define
            # the pattern-tree completely. We compute a number
            # that identifies this tree, and add 1 at the current pattern
            # in that tree's row.
            tree_id = get_pt_id(tree, v_idxs, perm_x, perm_y)
            if tree_id not in res:
                res[tree_id] = {}
            res[tree_id][p] = res[tree_id].get(p, 0) + 1

    return res.values()


def get_v_idxs(tree: Graph, comp: Composition) -> Dict[int, List[int]]:
    ''' Return a map of vertices to tree point indices. '''

    # The points are arbitrarily numbered 0,1,...,k-1.
    # No loss of generality since we iterate over all x-orderings above.
    res: Dict[int, List[int]] = {}
    idx = 0
    for v in range(len(tree)):
        for _ in range(comp[v]):
            if v not in res:
                res[v] = []
            res[v] += [idx]
            idx += 1
    return res


def allow_perm_x(v_idxs: Dict[int, List[int]], perm_x: Permutation) -> bool:
    ''' Check whether the x coordinates within each vertex are increasing. '''
    for idxs in v_idxs.values():
        for i in range(1, len(idxs)):
            idx = idxs[i]
            prev_idx = idxs[i - 1]
            if perm_x[idx] < perm_x[prev_idx]:
                return False
    return True


def get_pt_id(tree: Graph, v_idxs: Dict[int, List[int]],
              perm_x: Permutation, perm_y: Permutation) -> int:
    ''' Get an integer that identifies a pattern-tree given the tree
        structure, and the relative x and y ordering of its points. '''

    # We insert a bit into the result for every x- and y-ordering.
    # Note that we treat k, the graph structure, and the vertex sizes
    # as fixed (see gen_rows).
    res = 0

    # To have consistent IDs, it is crucial to iterate over the edges
    # in the same order every time.
    for u, w in tree.edge_iterator(labels=False, sort_vertices=True):
        assert u < w
        # In the edge {u, w}, every pair of pattern points from
        # u and v are relatively ordered in both x and y.
        # We add both bits to res.
        for u_idx in v_idxs[u]:
            for w_idx in v_idxs[w]:
                for p in perm_x, perm_y:
                    res <<= 1
                    res += int(p[u_idx] < p[w_idx])

    # Every pair of points inside a vertex are relatively ordered in
    # the y coordinate (the x coordinate is always increasing).
    for v in range(len(tree)):
        for idx1, idx2 in combinations(v_idxs[v], r=2):
            res <<= 1
            res += int(perm_y[idx1] < perm_y[idx2])
    return res


def gen_special5_rows() -> Iterable[PermInt]:
    ''' Generate rows according to the special 3214 and 43215 patterns,
        and their symmetries. '''

    # Insert the unweighted 3214 and 43215 patterns
    # that we can compute directly (see Section 4.3-4.5).
    row3214: PermInt = {Permutation([3, 2, 1, 4]): 1}
    row43215: PermInt = {Permutation([4, 3, 2, 1, 5]): 1}
    all_rows = [row3214, row43215]

    # Below we consider trees with weighted 3214 (see Section 4.2).
    # Since we are interested in 5-patterns, we only weigh a single point
    # in 3214 with a count of pairs. For example, consider the following tree:
    #
    #        ____________
    #       | root: 3214 |
    #        ‾‾‾‾‾‾‾‾‾‾‾‾
    #              |
    #              |  p1 = root2
    #              |
    #           ______
    #          | p: 1 |
    #           ‾‾‾‾‾‾
    #              |
    #              |  q1.x < p1.x
    #              |  q1.y > p1.y
    #              |
    #           ______
    #          | q: 1 |
    #           ‾‾‾‾‾‾
    #
    # The labels in this example signify that the root's weighted point is 2,
    # and the leaf-point q1 is in its top-left. We could have chosen any other
    # root-point and any other quadrant for q1's position relative to it,
    # corresponding to different trees.
    #
    # The linear combination of patterns counted by this tree can be computed
    # directly by inserting a point into 3214 in every possible position
    # in 2's top-left. For example, inserting the leaf-point left of 3 and
    # above 4 gives the pattern #53214.
    #
    # For the sake of code simplicity, we compute the rows in a different order
    # by iterating over patterns instead (similarly to gen_rows above).
    # Suppose we fix a 5-pattern, and in it fix a root-point and a leaf-point.
    # Suppose that omitting the leaf-point from the 5-pattern correctly
    # gives 3214.
    #
    # Given the root-point, the position of the leaf-point relative to it
    # fully determines the tree's labels. Hence we consider these the tree's
    # "ID" as above, and add 1 to its linear combination for this pattern.
    #
    # Note: a pair of unconstrained points are allowed to be equal,
    # so the above tree may also contain 3214 in the linear combination.
    # Since 3214 is spanned (Section 4.3), it can be subtracted from these
    # combinations, so we omit this case for simplicity.
    rows_dict = {}
    for perm5 in Permutations(5):
        for leaf_point in range(5):
            # The root permutation is obtained by omitting the leaf-point.
            perm4 = perm5[:leaf_point] + perm5[leaf_point + 1:]
            # Must be 3214 at the root.
            if not (perm4[2] < perm4[1] < perm4[0] < perm4[3]):
                continue
            for root_point in range(5):
                if root_point == leaf_point:
                    continue
                # Given a 5-pattern and the position of the root-point,
                # the pattern-tree is identified by its relative ordering
                # to the leaf-point.
                tree_id = (
                    root_point,
                    (root_point < leaf_point),                  # x-ordering
                    (perm5[root_point] < perm5[leaf_point]),    # y-ordering
                )
                if tree_id not in rows_dict:
                    rows_dict[tree_id] = {}
                cur_row = rows_dict[tree_id]
                cur_row[perm5] = cur_row.get(perm5, 0) + 1

    all_rows.extend(rows_dict.values())

    # For every row above, yield every symmetry
    # (by acting with D4 on the input permutation, see Section 2.1).
    for row in all_rows:
        for new_row in get_row_D4_orbit(row):
            yield new_row


def get_D4():
    ''' Get a list of Permutation functions that correspond to the
        8 symmetries of rotations and reflections.'''
    rev = Permutation.reverse
    inv = Permutation.inverse
    compl = Permutation.complement
    return [
        (lambda p: p), rev, compl, compose(rev, compl),
        inv, compose(rev, inv), compose(compl, inv),
        compose(rev, compose(compl, inv)),
    ]


def get_row_D4_orbit(row: PermInt) -> Iterable[PermInt]:
    ''' Given a sparse row that maps permutations to integers,
        return a sequence of rows corresponding to its 8 symmetries.'''
    for function in get_D4():
        new_row: PermInt = {}
        for p, val in row.items():
            q = function(p)
            new_row[q] = new_row.get(q, 0) + val
        yield new_row


def enum_perms(max_k: int) -> PermInt:
    ''' Return a mapping of S_{<=k} to successive integers. '''
    res: PermInt = {}
    for k in range(1, max_k + 1):
        for p in Permutations(k):
            res[p] = len(res)
    return res


if __name__ == '__main__':
    main()
