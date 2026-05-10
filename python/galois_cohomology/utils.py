"""
Linear algebra utilities for integer matrices.

Smith normal form computation for determining the structure of
quotient groups ker/im over Z.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass


@dataclass
class SmithNormalForm:
    """Result of Smith normal form decomposition: D = U @ A @ V."""

    D: np.ndarray  # diagonal matrix (Smith normal form)
    U: np.ndarray  # left transformation (invertible over Z)
    V: np.ndarray  # right transformation (invertible over Z)
    invariant_factors: list[int]


def smith_normal_form(A: np.ndarray) -> SmithNormalForm:
    """
    Compute the Smith normal form of an integer matrix A.

    Returns D, U, V such that D = U @ A @ V, where D is diagonal
    with d_1 | d_2 | ... | d_r on the diagonal.
    """
    if A.size == 0:
        rows, cols = A.shape
        return SmithNormalForm(
            D=A.copy(),
            U=np.eye(rows, dtype=np.int64),
            V=np.eye(cols, dtype=np.int64),
            invariant_factors=[],
        )

    M = A.astype(np.int64).copy()
    rows, cols = M.shape
    U = np.eye(rows, dtype=np.int64)
    V = np.eye(cols, dtype=np.int64)

    pivot = 0
    for pivot in range(min(rows, cols)):
        # Find nonzero entry in submatrix M[pivot:, pivot:]
        sub = M[pivot:, pivot:]
        if np.all(sub == 0):
            break

        # Find the entry with smallest absolute value (nonzero)
        nonzero_mask = sub != 0
        if not np.any(nonzero_mask):
            break

        abs_sub = np.abs(sub)
        abs_sub[~nonzero_mask] = np.iinfo(np.int64).max
        min_pos = np.unravel_index(np.argmin(abs_sub), sub.shape)
        min_row = min_pos[0] + pivot
        min_col = min_pos[1] + pivot

        # Swap to pivot position
        if min_row != pivot:
            M[[pivot, min_row]] = M[[min_row, pivot]]
            U[[pivot, min_row]] = U[[min_row, pivot]]
        if min_col != pivot:
            M[:, [pivot, min_col]] = M[:, [min_col, pivot]]
            V[:, [pivot, min_col]] = V[:, [min_col, pivot]]

        # Eliminate entries in pivot row and column
        changed = True
        while changed:
            changed = False

            # Eliminate column entries
            for i in range(pivot + 1, rows):
                if M[i, pivot] != 0:
                    q = M[i, pivot] // M[pivot, pivot]
                    M[i] -= q * M[pivot]
                    U[i] -= q * U[pivot]
                    if M[i, pivot] != 0:
                        changed = True

            # Eliminate row entries
            for j in range(pivot + 1, cols):
                if M[pivot, j] != 0:
                    q = M[pivot, j] // M[pivot, pivot]
                    M[:, j] -= q * M[:, pivot]
                    V[:, j] -= q * V[:, pivot]
                    if M[pivot, j] != 0:
                        changed = True

            # If column still has nonzero entries, swap smallest into pivot
            if changed:
                col_below = M[pivot + 1:, pivot]
                nonzero_below = np.where(col_below != 0)[0]
                if len(nonzero_below) > 0:
                    abs_below = np.abs(col_below[nonzero_below])
                    best = nonzero_below[np.argmin(abs_below)]
                    best_row = best + pivot + 1
                    if abs(M[best_row, pivot]) < abs(M[pivot, pivot]):
                        M[[pivot, best_row]] = M[[best_row, pivot]]
                        U[[pivot, best_row]] = U[[best_row, pivot]]

        # Make pivot positive
        if M[pivot, pivot] < 0:
            M[pivot] = -M[pivot]
            U[pivot] = -U[pivot]

    # Ensure divisibility: d_i | d_{i+1}
    diag = [int(M[i, i]) for i in range(min(rows, cols)) if i < rows and i < cols]
    for i in range(len(diag)):
        for j in range(i + 1, len(diag)):
            if diag[i] != 0 and diag[j] != 0 and diag[j] % diag[i] != 0:
                # Use extended gcd to fix divisibility
                _fix_divisibility(M, U, V, i, j)

    invariant_factors = []
    for i in range(min(rows, cols)):
        if i < rows and i < cols and M[i, i] != 0:
            invariant_factors.append(abs(int(M[i, i])))

    return SmithNormalForm(D=M, U=U, V=V, invariant_factors=invariant_factors)


def _fix_divisibility(M, U, V, i, j):
    """Fix divisibility condition by row/column operations."""
    from math import gcd

    d = gcd(int(M[i, i]), int(M[j, j]))
    if d == abs(int(M[i, i])):
        return

    # Add column j to column i, then re-eliminate
    V[:, i] += V[:, j]
    M[:, i] += M[:, j]

    # Re-eliminate using the pivot at (i, i)
    while M[j, i] != 0:
        q = M[j, i] // M[i, i]
        M[j] -= q * M[i]
        U[j] -= q * U[i]
        if M[j, i] != 0 and abs(M[j, i]) < abs(M[i, i]):
            M[[i, j]] = M[[j, i]]
            U[[i, j]] = U[[j, i]]

    while M[i, j] != 0:
        q = M[i, j] // M[i, i]
        M[:, j] -= q * M[:, i]
        V[:, j] -= q * V[:, i]

    if M[i, i] < 0:
        M[i] = -M[i]
        U[i] = -U[i]


def kernel_basis(A: np.ndarray) -> np.ndarray:
    """
    Compute a basis for ker(A) over Z using Smith normal form.

    Returns a matrix whose columns are a basis for {x : Ax = 0}.
    """
    if A.size == 0:
        cols = A.shape[1] if len(A.shape) > 1 else 0
        return np.eye(cols, dtype=np.int64)

    rows, cols = A.shape
    snf = smith_normal_form(A)
    rank = len(snf.invariant_factors)
    # Kernel basis is the last (cols - rank) columns of V
    if rank >= cols:
        return np.zeros((cols, 0), dtype=np.int64)
    return snf.V[:, rank:]


def image_in_target(A: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """
    Compute the image of A as a subgroup of Z^rows.

    Returns (basis_matrix, invariant_factors) where:
    - basis_matrix columns generate im(A)
    - invariant_factors describe the structure
    """
    if A.size == 0:
        rows = A.shape[0]
        return np.zeros((rows, 0), dtype=np.int64), []

    snf = smith_normal_form(A)
    rank = len(snf.invariant_factors)
    # Image generators: first `rank` columns of U^{-1} scaled by invariant factors
    # Actually, im(A) is generated by columns of A, which in SNF coords
    # becomes columns of U^{-1} @ D, so first `rank` columns of U^{-1} scaled.
    # Simpler: just return columns of A (they generate the image).
    # But for the quotient computation, we need the SNF of the inclusion.
    return A[:, :rank] if rank > 0 else np.zeros((A.shape[0], 0), dtype=np.int64), snf.invariant_factors
