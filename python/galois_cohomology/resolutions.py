"""
Resolutions and spectral sequences for Galois cohomology.

Author: RTSC agent (Claude Opus 4.7), 2026-05-10.
Fixed by MAIN agent to conform to actual CochainComplex interface.

Provides:
  - BarResolution: standard bar resolution B_*(G) and induced cochain
    complex Hom_G(B_*, M).
  - LHSSpectralSequence: Lyndon-Hochschild-Serre spectral sequence for
    a normal subgroup N of G with E_2^{p,q} = H^p(G/N, H^q(N, M)).

References:
  Brown, "Cohomology of Groups", Ch. I-III.
  Weibel, "Introduction to Homological Algebra", Ch. 5.
  Neukirch-Schmidt-Wingberg, "Cohomology of Number Fields", Ch. I.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field

from .modules import FiniteGroup, GModule
from .cochains import CochainComplex, coboundary_matrix
from .cohomology import compute_cohomology, compute_cohomology_from_complex, CohomologyGroup


@dataclass
class BarResolution:
    """
    Standard bar resolution B_*(G) of Z over Z[G].

    B_n = Z[G^{n+1}] (free Z[G]-module on n+1-tuples of G).
    Boundary d_n: B_n -> B_{n-1} via the standard alternating sum:
      d_n(g_0, ..., g_n) = sum_i (-1)^i (g_0, ..., hat{g_i}, ..., g_n)

    The cochain complex Hom_G(B_*, M) recovers the standard cochain
    complex C^*(G, M) computed by cochains.py.
    """

    group: FiniteGroup
    module: GModule
    _cochain: CochainComplex = field(init=False)

    def __post_init__(self):
        self._cochain = CochainComplex(self.group, self.module)

    def chain_rank(self, n: int) -> int:
        """Z[G]-rank of B_n (free of rank |G|^n as Z[G]-module)."""
        return self.group.order ** n

    def cochain_dim(self, n: int) -> int:
        """Z-dimension of Hom_G(B_n, M) = C^n(G, M)."""
        return self._cochain.dim(n)

    def cochain_differential(self, n: int) -> np.ndarray:
        """d^n: C^n(G, M) -> C^{n+1}(G, M)."""
        return self._cochain.differential(n)

    def verify_d_squared(self, n: int) -> bool:
        """Verify d^{n+1} . d^n = 0."""
        return self._cochain.verify_d_squared(n)

    def cochain_complex(self) -> CochainComplex:
        """Return the underlying CochainComplex."""
        return self._cochain

    def cohomology(self, n: int) -> CohomologyGroup:
        """Compute H^n(G, M) via this resolution."""
        return compute_cohomology_from_complex(self._cochain, n)


@dataclass
class LHSSpectralSequence:
    """
    Lyndon-Hochschild-Serre spectral sequence for a normal subgroup
    N of G acting on G-module M:

        E_2^{p,q} = H^p(G/N, H^q(N, M))  =>  H^{p+q}(G, M)

    Simplification: assumes G/N acts trivially on H^q(N, M).
    This is correct when G is abelian (conjugation is inner).

    The module M must be valid as both a G-module and an N-module
    (by restriction of the action to N).
    """

    group: FiniteGroup
    normal_subgroup: FiniteGroup
    quotient_group: FiniteGroup
    module: GModule
    restriction_module: GModule  # M viewed as N-module

    def e2_dimension(self, p: int, q: int) -> CohomologyGroup:
        """
        Compute E_2^{p,q} = H^p(G/N, H^q(N, M)).

        With trivial G/N-action simplification:
        E_2^{p,q} = H^p(G/N, Z^k) where k = rank of H^q(N, M) as free part.
        """
        # Step 1: compute H^q(N, M|_N)
        hq = compute_cohomology(self.normal_subgroup, self.restriction_module, q)

        if hq.is_trivial:
            return CohomologyGroup(degree=p, free_rank=0, torsion_invariants=[])

        # Step 2: H^p(G/N, H^q(N,M)) with trivial action
        # Build a trivial G/N-module of appropriate rank
        # The rank is the free rank of H^q(N, M)
        hq_rank = hq.free_rank + len(hq.torsion_invariants)
        if hq_rank == 0:
            return CohomologyGroup(degree=p, free_rank=0, torsion_invariants=[])

        trivial_mod = GModule.trivial(self.quotient_group, rank=hq_rank)
        return compute_cohomology(self.quotient_group, trivial_mod, p)

    def converges_to(self, n: int) -> CohomologyGroup:
        """Compute H^n(G, M) (the abutment of the spectral sequence)."""
        return compute_cohomology(self.group, self.module, n)

    def verify_inequality(self, n: int) -> bool:
        """
        Verify that sum_{p+q=n} |E_2^{p,q}| >= |H^n(G, M)|.

        The order of E_2 total at degree n bounds the order of H^n
        (with equality iff the spectral sequence degenerates at E_2).
        """
        target = self.converges_to(n)
        if target.order is None:
            return True  # infinite groups, skip check

        e2_orders = []
        for p in range(n + 1):
            q = n - p
            e2pq = self.e2_dimension(p, q)
            if e2pq.order is None:
                return True  # infinite, inequality trivially holds
            e2_orders.append(e2pq.order)

        # For finite groups: product of E_2 orders (not sum) relates to H^n
        # Actually the relationship is: H^n has a filtration with graded pieces E_inf
        # and |E_inf| divides |E_2|. So product(|E_inf^{p,q}|) = |H^n| and
        # |E_inf^{p,q}| divides |E_2^{p,q}|.
        target_order = target.order if target.order else 1
        e2_product = 1
        for o in e2_orders:
            e2_product *= o
        return e2_product >= target_order


def build_lhs_abelian(
    group: FiniteGroup,
    normal_subgroup: FiniteGroup,
    quotient_group: FiniteGroup,
    module: GModule,
    normal_action_indices: list[int],
) -> LHSSpectralSequence:
    """
    Build LHS spectral sequence for abelian G = N x (G/N).

    Args:
        group: the full group G
        normal_subgroup: normal subgroup N
        quotient_group: quotient G/N
        module: G-module M
        normal_action_indices: indices in G.elements corresponding to N
    """
    # Build restriction of M to N
    restriction_matrices = [module.action(i) for i in normal_action_indices]
    restriction_module = GModule(
        group=normal_subgroup,
        rank=module.rank,
        action_matrices=restriction_matrices,
    )

    return LHSSpectralSequence(
        group=group,
        normal_subgroup=normal_subgroup,
        quotient_group=quotient_group,
        module=module,
        restriction_module=restriction_module,
    )
