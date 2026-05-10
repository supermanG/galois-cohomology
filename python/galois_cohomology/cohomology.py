"""
Cohomology group computation: H^n(G, M) = ker(d^n) / im(d^{n-1}).

For finite groups G acting on finitely generated Z-modules M,
H^n(G, M) is a finitely generated abelian group. We compute its
structure (free rank + torsion invariant factors) via Smith normal form.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from math import gcd
from functools import reduce

from .modules import FiniteGroup, GModule
from .cochains import CochainComplex, coboundary_matrix
from .utils import smith_normal_form, kernel_basis


@dataclass
class CohomologyGroup:
    """
    Structure of H^n(G, M) as a finitely generated abelian group.

    H^n = Z^free_rank + Z/d_1 + Z/d_2 + ... + Z/d_k
    where d_1 | d_2 | ... | d_k.
    """

    degree: int
    free_rank: int
    torsion_invariants: list[int]  # d_i > 1, with d_1 | d_2 | ...
    cocycle_basis: np.ndarray | None = None  # columns span ker(d^n)
    coboundary_basis: np.ndarray | None = None  # columns span im(d^{n-1})

    @property
    def order(self) -> int | None:
        """Order of the group if finite, None if infinite."""
        if self.free_rank > 0:
            return None
        if not self.torsion_invariants:
            return 1
        result = 1
        for d in self.torsion_invariants:
            result *= d
        return result

    @property
    def is_trivial(self) -> bool:
        return self.free_rank == 0 and len(self.torsion_invariants) == 0

    def __repr__(self) -> str:
        parts = []
        if self.free_rank > 0:
            parts.append(f"Z^{self.free_rank}" if self.free_rank > 1 else "Z")
        for d in self.torsion_invariants:
            parts.append(f"Z/{d}Z")
        if not parts:
            return f"H^{self.degree} = 0"
        return f"H^{self.degree} = {' + '.join(parts)}"


def compute_cohomology(
    group: FiniteGroup, module: GModule, n: int, keep_bases: bool = False
) -> CohomologyGroup:
    """
    Compute H^n(G, M) = ker(d^n) / im(d^{n-1}).

    Args:
        group: finite group G
        module: G-module M
        n: cohomological degree
        keep_bases: if True, store cocycle/coboundary bases in result

    Returns:
        CohomologyGroup describing the structure of H^n(G, M)
    """
    complex_ = CochainComplex(group, module)
    return compute_cohomology_from_complex(complex_, n, keep_bases)


def compute_cohomology_from_complex(
    complex_: CochainComplex, n: int, keep_bases: bool = False
) -> CohomologyGroup:
    """Compute H^n from a cochain complex."""
    dn = complex_.differential(n)
    cocycles = kernel_basis(dn)  # columns span Z^n = ker(d^n)

    if n == 0:
        # im(d^{-1}) = 0
        coboundaries = np.zeros((complex_.dim(0), 0), dtype=np.int64)
    else:
        dn_minus_1 = complex_.differential(n - 1)
        coboundaries = dn_minus_1  # columns of d^{n-1} span im(d^{n-1}) in C^n

    # Compute H^n = ker(d^n) / im(d^{n-1})
    # Express coboundaries in the cocycle basis, then take quotient
    free_rank, torsion = _quotient_structure(cocycles, coboundaries)

    return CohomologyGroup(
        degree=n,
        free_rank=free_rank,
        torsion_invariants=torsion,
        cocycle_basis=cocycles if keep_bases else None,
        coboundary_basis=coboundaries if keep_bases else None,
    )


def _quotient_structure(
    ker_basis: np.ndarray, im_generators: np.ndarray
) -> tuple[int, list[int]]:
    """
    Compute the structure of ker/im as a finitely generated abelian group.

    ker_basis: columns form a Z-basis of ker (a free Z-module)
    im_generators: columns generate im as a subgroup of the ambient Z^n

    Strategy:
    1. Express im generators in the ker basis coordinates.
    2. Take Smith normal form of the coordinate matrix.
    3. Read off the quotient structure.
    """
    ker_rank = ker_basis.shape[1] if ker_basis.ndim == 2 else 0

    if ker_rank == 0:
        return 0, []

    if im_generators.size == 0 or im_generators.shape[1] == 0:
        return ker_rank, []

    # Express im generators in ker_basis coordinates.
    # If ker_basis is K (ambient_dim x ker_rank), we need to solve K @ x = b
    # for each column b of im_generators that lies in the column space of K.
    # Since K has full column rank (it's a basis), we use least squares over Z.
    coord_matrix = _express_in_basis(ker_basis, im_generators)

    if coord_matrix is None or coord_matrix.size == 0:
        return ker_rank, []

    # Smith normal form of the coordinate matrix gives the quotient structure
    snf = smith_normal_form(coord_matrix)
    torsion = [d for d in snf.invariant_factors if d > 1]
    rank_of_im_in_ker = len(snf.invariant_factors)
    free_rank = ker_rank - rank_of_im_in_ker

    return free_rank, torsion


def _express_in_basis(basis: np.ndarray, vectors: np.ndarray) -> np.ndarray | None:
    """
    Express vectors as integer linear combinations of basis columns.

    Returns coordinate matrix C such that basis @ C = vectors (approximately,
    for those vectors that lie in the span).
    """
    if basis.shape[1] == 0:
        return None

    ambient_dim, ker_rank = basis.shape
    _, num_vectors = vectors.shape

    # Use least-squares and round (works when vectors are truly in the span)
    # For exact integer computation, we solve via the normal equations
    # B^T B x = B^T v for each column v
    BtB = basis.T @ basis
    Btv = basis.T @ vectors

    # Solve using sympy for exact integer arithmetic
    try:
        from sympy import Matrix

        BtB_sym = Matrix(BtB.tolist())
        Btv_sym = Matrix(Btv.tolist())

        result = np.zeros((ker_rank, num_vectors), dtype=np.int64)
        for j in range(num_vectors):
            col = Btv_sym[:, j]
            sol = BtB_sym.solve(col)
            for i in range(ker_rank):
                result[i, j] = int(sol[i])
        return result
    except Exception:
        # Fallback: numpy least squares with rounding
        result, _, _, _ = np.linalg.lstsq(basis.astype(float), vectors.astype(float), rcond=None)
        return np.round(result).astype(np.int64)


def fixed_points(group: FiniteGroup, module: GModule) -> np.ndarray:
    """
    Compute M^G = H^0(G, M), the submodule of G-fixed elements.

    Returns a matrix whose columns are a basis for M^G.
    """
    # M^G = {m in M : g.m = m for all g in G}
    # This is the intersection of ker(action(g) - I) for all g
    r = module.rank
    constraints = []
    for g in range(group.order):
        diff = module.action(g) - np.eye(r, dtype=np.int64)
        constraints.append(diff)

    if not constraints:
        return np.eye(r, dtype=np.int64)

    # Stack all constraints
    big_matrix = np.vstack(constraints)
    return kernel_basis(big_matrix)
