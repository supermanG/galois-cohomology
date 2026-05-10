"""
Cochain complexes and coboundary operators for group cohomology.

C^n(G, M) = Map(G^n, M) with the standard coboundary:

(d^n f)(g_0, ..., g_n) = g_0 . f(g_1, ..., g_n)
    + sum_{i=1}^{n} (-1)^i f(g_0, ..., g_{i-1}*g_i, ..., g_n)
    + (-1)^{n+1} f(g_0, ..., g_{n-1})

For finite G with |G| = m and M = Z^r, we represent C^n(G, M) as Z^{r * m^n}
and d^n as an integer matrix.
"""

from __future__ import annotations

import numpy as np
from itertools import product as cartesian_product
from dataclasses import dataclass
from typing import Optional

from .modules import FiniteGroup, GModule


def _multi_index(indices: tuple[int, ...], base: int) -> int:
    """Convert a tuple of group element indices to a flat index in G^n."""
    result = 0
    for idx in indices:
        result = result * base + idx
    return result


def _from_multi_index(flat: int, n: int, base: int) -> tuple[int, ...]:
    """Convert flat index back to tuple of group element indices."""
    indices = []
    for _ in range(n):
        indices.append(flat % base)
        flat //= base
    return tuple(reversed(indices))


@dataclass
class CochainGroup:
    """
    C^n(G, M): the group of n-cochains.

    Represented as Z^dim where dim = rank(M) * |G|^n.
    Basis ordering: for f in C^n(G, M), the coefficient of the basis vector
    e_j at position (g_0, ..., g_{n-1}) is at flat index:
        _multi_index((g_0, ..., g_{n-1}), |G|) * rank + j
    """

    group: FiniteGroup
    module: GModule
    degree: int

    @property
    def dim(self) -> int:
        return self.module.rank * (self.group.order ** self.degree)


def coboundary_matrix(group: FiniteGroup, module: GModule, n: int) -> np.ndarray:
    """
    Compute the matrix of d^n: C^n(G, M) -> C^{n+1}(G, M).

    The coboundary formula (inhomogeneous bar notation):
    (d^n f)(g_0, ..., g_n) = g_0 . f(g_1, ..., g_n)
        + sum_{i=1}^{n} (-1)^i f(g_0, ..., g_{i-1}*g_i, ..., g_n)
        + (-1)^{n+1} f(g_0, ..., g_{n-1})
    """
    m = group.order
    r = module.rank
    source_dim = r * (m ** n)
    target_dim = r * (m ** (n + 1))

    if n == 0:
        # d^0: C^0(G, M) -> C^1(G, M)
        # (d^0 f)(g) = g.f() - f()
        # C^0 = M = Z^r, C^1 = Map(G, M) = Z^{r*m}
        mat = np.zeros((target_dim, source_dim), dtype=np.int64)
        for g in range(m):
            row_offset = g * r
            action_g = module.action(g)
            for j in range(r):
                for k in range(r):
                    mat[row_offset + j, k] += int(action_g[j, k])
                mat[row_offset + j, j] -= 1
        return mat

    mat = np.zeros((target_dim, source_dim), dtype=np.int64)

    for target_tuple in cartesian_product(range(m), repeat=n + 1):
        g = target_tuple
        target_base = _multi_index(g, m) * r

        # Term 0: g_0 . f(g_1, ..., g_n)
        source_tuple_0 = g[1:]
        source_base_0 = _multi_index(source_tuple_0, m) * r
        action_g0 = module.action(g[0])
        for j in range(r):
            for k in range(r):
                mat[target_base + j, source_base_0 + k] += int(action_g0[j, k])

        # Terms i=1..n: (-1)^i f(g_0, ..., g_{i-1}*g_i, ..., g_n)
        for i in range(1, n + 1):
            sign = (-1) ** i
            merged = group.mult(g[i - 1], g[i])
            source_tuple_i = g[:i-1] + (merged,) + g[i+1:]
            source_base_i = _multi_index(source_tuple_i, m) * r
            for j in range(r):
                mat[target_base + j, source_base_i + j] += sign

        # Term n+1: (-1)^{n+1} f(g_0, ..., g_{n-1})
        sign = (-1) ** (n + 1)
        source_tuple_last = g[:n]
        source_base_last = _multi_index(source_tuple_last, m) * r
        for j in range(r):
            mat[target_base + j, source_base_last + j] += sign

    return mat


@dataclass
class CochainComplex:
    """
    The standard cochain complex C^*(G, M).

    Lazily computes and caches coboundary matrices.
    """

    group: FiniteGroup
    module: GModule
    _differentials: dict[int, np.ndarray]

    def __init__(self, group: FiniteGroup, module: GModule):
        self.group = group
        self.module = module
        self._differentials = {}

    def cochain_group(self, n: int) -> CochainGroup:
        return CochainGroup(group=self.group, module=self.module, degree=n)

    def differential(self, n: int) -> np.ndarray:
        """Return d^n: C^n -> C^{n+1}."""
        if n < 0:
            target_dim = self.module.rank * (self.group.order ** (n + 1))
            source_dim = 0
            return np.zeros((target_dim, source_dim), dtype=np.int64)
        if n not in self._differentials:
            self._differentials[n] = coboundary_matrix(self.group, self.module, n)
        return self._differentials[n]

    def verify_d_squared(self, n: int) -> bool:
        """Verify d^{n+1} . d^n = 0."""
        dn = self.differential(n)
        dn1 = self.differential(n + 1)
        product = dn1 @ dn
        return np.all(product == 0)

    def dim(self, n: int) -> int:
        """Dimension of C^n(G, M)."""
        if n < 0:
            return 0
        return self.module.rank * (self.group.order ** n)
