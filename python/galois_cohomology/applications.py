"""
Classical applications of Galois cohomology.

- Hilbert's Theorem 90: H^1(Gal(L/K), L*) = 0
- Fixed point computation for Galois descent
- Norm maps and Herbrand quotient setup
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from .modules import FiniteGroup, GModule
from .cohomology import compute_cohomology, CohomologyGroup, fixed_points
from .fields import QuadraticField, CyclotomicField


def verify_hilbert_90_quadratic(d: int) -> CohomologyGroup:
    """
    Verify Hilbert 90 for Q(sqrt(d))/Q.

    H^1(C_2, M) where M is the additive Galois module Q(sqrt(d)).
    For the additive version, H^1(Gal, L) = 0 (additive Hilbert 90).

    Note: The multiplicative Hilbert 90 (H^1(Gal, L*) = 0) requires
    the multiplicative group which is not finitely generated over Z.
    We verify the additive analogue.
    """
    K = QuadraticField(d=d)
    G = K.galois_group()
    M = K.multiplicative_module()  # rank-2 Z-module with C_2 action
    return compute_cohomology(G, M, 1)


def verify_hilbert_90_cyclic(n: int) -> CohomologyGroup:
    """
    Verify additive Hilbert 90 for Q(zeta_n)/Q.

    H^1(Gal(Q(zeta_n)/Q), Z^{phi(n)}) where the action is the Galois action.
    """
    F = CyclotomicField(n)
    G = F.galois_group()
    M = F.standard_module()
    return compute_cohomology(G, M, 1)


def compute_h2_galois(n: int) -> CohomologyGroup:
    """
    Compute H^2(Gal(Q(zeta_n)/Q), standard module).

    This gives part of the Brauer group computation.
    """
    F = CyclotomicField(n)
    G = F.galois_group()
    M = F.standard_module()
    return compute_cohomology(G, M, 2)


def galois_descent_quadratic(d: int) -> np.ndarray:
    """
    Compute the descended subspace (Q(sqrt(d)))^{Gal} = Q.

    Returns basis for the fixed submodule M^G.
    """
    K = QuadraticField(d=d)
    G = K.galois_group()
    M = K.multiplicative_module()
    return fixed_points(G, M)


@dataclass
class NormMap:
    """
    The norm map N: M -> M^G defined by N(m) = sum_{g in G} g.m.

    For Tate cohomology, we need both the norm and its kernel.
    """

    group: FiniteGroup
    module: GModule
    matrix: np.ndarray

    @classmethod
    def build(cls, group: FiniteGroup, module: GModule) -> NormMap:
        """Construct the norm map N = sum_g action(g)."""
        r = module.rank
        norm_mat = np.zeros((r, r), dtype=np.int64)
        for g in range(group.order):
            norm_mat += module.action(g)
        return cls(group=group, module=module, matrix=norm_mat)

    def kernel_rank(self) -> int:
        """Rank of ker(N)."""
        from .utils import kernel_basis
        ker = kernel_basis(self.matrix)
        return ker.shape[1] if ker.ndim == 2 else 0

    def image_rank(self) -> int:
        """Rank of im(N)."""
        from .utils import smith_normal_form
        snf = smith_normal_form(self.matrix)
        return len(snf.invariant_factors)


def norm_map(group: FiniteGroup, module: GModule) -> NormMap:
    """Construct the norm map for the given G-module."""
    return NormMap.build(group, module)
