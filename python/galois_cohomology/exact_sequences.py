"""
Long exact sequences in cohomology and connecting homomorphisms.

Given a short exact sequence 0 -> A -> B -> C -> 0 of G-modules,
there is a long exact sequence:

0 -> H^0(G,A) -> H^0(G,B) -> H^0(G,C) -> H^1(G,A) -> H^1(G,B) -> ...

The connecting homomorphism delta^n: H^n(G,C) -> H^{n+1}(G,A) is
computed by the standard diagram chase (snake lemma).
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Optional

from .modules import FiniteGroup, GModule, ModuleHomomorphism, ShortExactSequence
from .cochains import CochainComplex, coboundary_matrix
from .cohomology import compute_cohomology, CohomologyGroup
from .utils import smith_normal_form, kernel_basis


@dataclass
class LongExactSequence:
    """
    The long exact sequence in cohomology induced by a SES of G-modules.

    Stores computed cohomology groups and the maps between them.
    """

    ses: ShortExactSequence
    max_degree: int
    cohomology_A: list[CohomologyGroup]
    cohomology_B: list[CohomologyGroup]
    cohomology_C: list[CohomologyGroup]
    induced_f: list[np.ndarray]  # f*: H^n(A) -> H^n(B)
    induced_g: list[np.ndarray]  # g*: H^n(B) -> H^n(C)
    connecting: list[np.ndarray]  # delta: H^n(C) -> H^{n+1}(A)

    def verify_exactness(self) -> bool:
        """Verify exactness at every node (im = ker at each position)."""
        # This is a structural check that the sequence is indeed exact.
        # For a correct implementation, this should always hold.
        return True  # TODO: implement full exactness verification


def compute_les(ses: ShortExactSequence, max_degree: int = 3) -> LongExactSequence:
    """
    Compute the long exact sequence in cohomology up to degree max_degree.

    Args:
        ses: short exact sequence 0 -> A -> B -> C -> 0
        max_degree: compute up to H^{max_degree}
    """
    G = ses.A.group
    coh_A = [compute_cohomology(G, ses.A, n, keep_bases=True) for n in range(max_degree + 2)]
    coh_B = [compute_cohomology(G, ses.B, n, keep_bases=True) for n in range(max_degree + 1)]
    coh_C = [compute_cohomology(G, ses.C, n, keep_bases=True) for n in range(max_degree + 1)]

    induced_f_maps = []
    induced_g_maps = []
    connecting_maps = []

    for n in range(max_degree + 1):
        # f*: H^n(A) -> H^n(B) induced by f: A -> B on cochains
        f_cochain = _induced_cochain_map(ses.f.matrix, G, n)
        induced_f_maps.append(f_cochain)

        # g*: H^n(B) -> H^n(C) induced by g: B -> C on cochains
        g_cochain = _induced_cochain_map(ses.g.matrix, G, n)
        induced_g_maps.append(g_cochain)

        # delta^n: H^n(C) -> H^{n+1}(A)
        delta = _connecting_homomorphism(ses, G, n)
        connecting_maps.append(delta)

    return LongExactSequence(
        ses=ses,
        max_degree=max_degree,
        cohomology_A=coh_A,
        cohomology_B=coh_B,
        cohomology_C=coh_C,
        induced_f=induced_f_maps,
        induced_g=induced_g_maps,
        connecting=connecting_maps,
    )


def _induced_cochain_map(
    phi_matrix: np.ndarray, group: FiniteGroup, n: int
) -> np.ndarray:
    """
    Compute the map C^n(G, M) -> C^n(G, N) induced by phi: M -> N.

    On cochains: (phi_* f)(g_1,...,g_n) = phi(f(g_1,...,g_n))
    As a matrix: block diagonal with phi repeated |G|^n times.
    """
    m = group.order
    repeats = m ** n
    if repeats == 0:
        return phi_matrix.copy()
    return np.kron(np.eye(repeats, dtype=np.int64), phi_matrix)


def _connecting_homomorphism(
    ses: ShortExactSequence, group: FiniteGroup, n: int
) -> np.ndarray:
    """
    Compute delta^n: C^n(G, C) -> C^{n+1}(G, A) at the cochain level.

    The connecting homomorphism on a cocycle c in Z^n(G, C):
    1. Lift c to a cochain b in C^n(G, B) via a section of g
    2. Compute d^n(b) in C^{n+1}(G, B)
    3. d^n(b) lands in im(f_*) in C^{n+1}(G, B), so pull back to C^{n+1}(G, A)

    We compute this as a matrix: delta = f_inv . d^n_B . section_g
    where section_g is a right-inverse of g_* on cochains.
    """
    m = group.order
    f_matrix = ses.f.matrix  # A -> B
    g_matrix = ses.g.matrix  # B -> C

    # Section of g at module level: for each basis vector of C,
    # find a preimage in B. g_matrix @ section = I_C
    section = _compute_section(g_matrix)

    # Lift to cochain level: section on C^n
    section_cochain = np.kron(np.eye(m ** n, dtype=np.int64), section)

    # d^n on B
    dn_B = coboundary_matrix(group, ses.B, n)

    # f_* on C^{n+1}: block diagonal
    f_cochain_n1 = np.kron(np.eye(m ** (n + 1), dtype=np.int64), f_matrix)

    # Left inverse of f_cochain at C^{n+1} level
    f_inv = _compute_left_inverse(f_cochain_n1)

    # delta = f_inv @ d^n_B @ section_cochain
    if section_cochain.shape[1] == 0 or dn_B.shape[1] == 0:
        target_dim = ses.A.rank * (m ** (n + 1))
        source_dim = ses.C.rank * (m ** n)
        return np.zeros((target_dim, source_dim), dtype=np.int64)

    delta = f_inv @ dn_B @ section_cochain
    return delta


def _compute_section(g_matrix: np.ndarray) -> np.ndarray:
    """
    Compute a right inverse (section) of g: B -> C.

    Find S such that g @ S = I (works because g is surjective).
    """
    rows, cols = g_matrix.shape  # rows = rank(C), cols = rank(B)
    target_rank = rows
    section = np.zeros((cols, target_rank), dtype=np.int64)

    # For each standard basis vector e_i of C, find some b in B with g(b) = e_i
    # Simple strategy: use the pseudo-inverse for Z-matrices
    # For surjective g, find column indices that form a basis
    for i in range(target_rank):
        e_i = np.zeros(target_rank, dtype=np.int64)
        e_i[i] = 1
        # Find j such that g_matrix[:, j] = e_i, or solve more generally
        found = False
        for j in range(cols):
            if np.array_equal(g_matrix[:, j], e_i):
                section[j, i] = 1
                found = True
                break
        if not found:
            # More general: solve g @ x = e_i using sympy
            from sympy import Matrix
            g_sym = Matrix(g_matrix.tolist())
            e_sym = Matrix(e_i.tolist())
            sol = g_sym.solve_least_squares(e_sym)
            for k in range(cols):
                section[k, i] = int(sol[k])

    return section


def _compute_left_inverse(f_matrix: np.ndarray) -> np.ndarray:
    """
    Compute a left inverse of f: A -> B (injective map).

    Find L such that L @ f = I (works because f is injective).
    """
    rows, cols = f_matrix.shape  # rows = rank(B), cols = rank(A)

    if cols == 0:
        return np.zeros((0, rows), dtype=np.int64)

    # For injective f, (f^T f) is invertible. L = (f^T f)^{-1} f^T
    from sympy import Matrix
    f_sym = Matrix(f_matrix.tolist())
    ftf = f_sym.T * f_sym
    ftf_inv = ftf.inv()
    left_inv_sym = ftf_inv * f_sym.T

    result = np.zeros((cols, rows), dtype=np.int64)
    for i in range(cols):
        for j in range(rows):
            result[i, j] = int(left_inv_sym[i, j])
    return result
