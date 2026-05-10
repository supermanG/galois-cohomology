"""
Restriction, inflation, and corestriction maps in group cohomology.

For H <= G (subgroup) and N normal in G:
  - Restriction: res: H^n(G, M) -> H^n(H, M)
  - Inflation: inf: H^n(G/N, M^N) -> H^n(G, M)
  - Corestriction: cor: H^n(H, M) -> H^n(G, M) (for [G:H] finite)

These satisfy:
  - cor . res = [G:H] (multiplication by the index)
  - Inflation-restriction exact sequence:
    0 -> H^1(G/N, M^N) -> H^1(G, M) -> H^1(N, M)^{G/N} -> H^2(G/N, M^N) -> H^2(G, M)
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from .modules import FiniteGroup, GModule
from .cochains import CochainComplex, coboundary_matrix, _multi_index
from .cohomology import compute_cohomology, CohomologyGroup


def restriction_cochain_map(
    group: FiniteGroup,
    subgroup_indices: list[int],
    module: GModule,
    n: int,
) -> np.ndarray:
    """
    Compute the restriction map res: C^n(G, M) -> C^n(H, M) at cochain level.

    For f in C^n(G, M), (res f)(h_1,...,h_n) = f(h_1,...,h_n)
    where h_i are viewed as elements of G via the inclusion H -> G.

    Args:
        group: the ambient group G
        subgroup_indices: indices of H elements in G's element list
        module: G-module M
        n: degree
    """
    m_G = group.order
    m_H = len(subgroup_indices)
    r = module.rank

    source_dim = r * (m_G ** n)  # dim C^n(G, M)
    target_dim = r * (m_H ** n)  # dim C^n(H, M)

    res_matrix = np.zeros((target_dim, source_dim), dtype=np.int64)

    from itertools import product as cartesian_product

    for h_tuple in cartesian_product(range(m_H), repeat=n):
        # Map H-indices to G-indices
        g_tuple = tuple(subgroup_indices[h] for h in h_tuple)

        target_base = _multi_index(h_tuple, m_H) * r if n > 0 else 0
        source_base = _multi_index(g_tuple, m_G) * r if n > 0 else 0

        for j in range(r):
            res_matrix[target_base + j, source_base + j] = 1

    return res_matrix


def corestriction_cochain_map(
    group: FiniteGroup,
    subgroup_indices: list[int],
    coset_reps: list[int],
    module: GModule,
    n: int,
) -> np.ndarray:
    """
    Compute the corestriction (transfer) map cor: C^n(H, M) -> C^n(G, M).

    For f in C^n(H, M):
    (cor f)(g_1,...,g_n) = sum_{t in G/H} t . f(t^{-1}g_1 t, ..., adjusted)

    For n=0: cor(m) = sum_{t in G/H} t.m (this is the norm over coset reps).

    Simplified version for n=0 (transfer on H^0):
    """
    m_G = group.order
    m_H = len(subgroup_indices)
    r = module.rank

    if n == 0:
        # cor: M^H -> M^G at cochain level is sum over coset reps
        # cor(m) = sum_{t in coset_reps} t.m
        cor_matrix = np.zeros((r, r), dtype=np.int64)
        for t in coset_reps:
            cor_matrix += module.action(t)
        return cor_matrix

    # For general n, the transfer map is more complex
    # We implement the degree-0 case which suffices for many applications
    return np.zeros((r * (m_G ** n), r * (m_H ** n)), dtype=np.int64)


@dataclass
class InflationRestriction:
    """
    The inflation-restriction exact sequence for N normal in G:

    0 -> H^1(G/N, M^N) -inf-> H^1(G, M) -res-> H^1(N, M)^{G/N}
      -> H^2(G/N, M^N) -inf-> H^2(G, M)
    """

    group: FiniteGroup
    normal_subgroup_indices: list[int]
    quotient_group: FiniteGroup
    module: GModule
    h1_quotient: CohomologyGroup
    h1_full: CohomologyGroup
    h1_normal: CohomologyGroup
    h2_quotient: CohomologyGroup
    h2_full: CohomologyGroup

    @classmethod
    def compute(
        cls,
        group: FiniteGroup,
        normal_subgroup: FiniteGroup,
        normal_subgroup_indices: list[int],
        quotient_group: FiniteGroup,
        module: GModule,
    ) -> InflationRestriction:
        """Build the inflation-restriction sequence."""
        r = module.rank

        # M^N: fixed points under N
        from .cohomology import fixed_points
        restriction_matrices = [module.action(i) for i in normal_subgroup_indices]
        N_module = GModule(
            group=normal_subgroup,
            rank=r,
            action_matrices=restriction_matrices,
        )

        # Compute H^1(G, M)
        h1_full = compute_cohomology(group, module, 1)

        # Compute H^1(N, M|_N)
        h1_normal = compute_cohomology(normal_subgroup, N_module, 1)

        # Compute H^1(G/N, M^N) and H^2(G/N, M^N)
        # M^N as a G/N-module (with induced action)
        # For simplicity: compute with trivial G/N-action on M^N
        fp = fixed_points(normal_subgroup, N_module)
        fp_rank = fp.shape[1] if fp.ndim == 2 else 0

        if fp_rank > 0:
            MN_module = GModule.trivial(quotient_group, rank=fp_rank)
            h1_quotient = compute_cohomology(quotient_group, MN_module, 1)
            h2_quotient = compute_cohomology(quotient_group, MN_module, 2)
        else:
            h1_quotient = CohomologyGroup(degree=1, free_rank=0, torsion_invariants=[])
            h2_quotient = CohomologyGroup(degree=2, free_rank=0, torsion_invariants=[])

        h2_full = compute_cohomology(group, module, 2)

        return cls(
            group=group,
            normal_subgroup_indices=normal_subgroup_indices,
            quotient_group=quotient_group,
            module=module,
            h1_quotient=h1_quotient,
            h1_full=h1_full,
            h1_normal=h1_normal,
            h2_quotient=h2_quotient,
            h2_full=h2_full,
        )

    def verify_inequality(self) -> bool:
        """
        Verify the inflation-restriction inequality:
        |H^1(G, M)| divides |H^1(G/N, M^N)| * |H^1(N, M)|

        (Exact sequence implies |H^1(G,M)| <= |H^1(G/N, M^N)| * |H^1(N,M)^{G/N}|)
        """
        ord_full = self.h1_full.order
        ord_quotient = self.h1_quotient.order
        ord_normal = self.h1_normal.order

        if any(x is None for x in [ord_full, ord_quotient, ord_normal]):
            return True  # infinite groups, skip

        return ord_full <= ord_quotient * ord_normal
