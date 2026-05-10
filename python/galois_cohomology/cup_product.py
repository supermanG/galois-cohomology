"""
Cup product on group cohomology.

The cup product H^p(G, M) x H^q(G, N) -> H^{p+q}(G, M tensor N)
is defined at the cochain level by:

  (f cup g)(g_1, ..., g_{p+q}) = f(g_1, ..., g_p) tensor g_1...g_p . g(g_{p+1}, ..., g_{p+q})

For M = N = Z with trivial action, this gives H^*(G, Z) a graded ring.
"""

from __future__ import annotations

import numpy as np
from itertools import product as cartesian_product

from .modules import FiniteGroup, GModule
from .cochains import CochainComplex, _multi_index, _from_multi_index


def cup_product_cochain(
    group: FiniteGroup,
    module_m: GModule,
    module_n: GModule,
    f: np.ndarray,
    g: np.ndarray,
    p: int,
    q: int,
) -> np.ndarray:
    """
    Compute the cup product f cup g at the cochain level.

    f in C^p(G, M), g in C^q(G, N) -> f cup g in C^{p+q}(G, M tensor N).

    (f cup g)(g_1,...,g_{p+q}) = f(g_1,...,g_p) tensor (g_1...g_p).g(g_{p+1},...,g_{p+q})

    For rank-1 trivial modules (M = N = Z), the tensor product is just multiplication.
    """
    m = group.order
    r_m = module_m.rank
    r_n = module_n.rank
    r_tensor = r_m * r_n  # rank of M tensor N

    result_dim = r_tensor * (m ** (p + q))
    result = np.zeros(result_dim, dtype=np.int64)

    for target_tuple in cartesian_product(range(m), repeat=p + q):
        g_all = target_tuple
        g_first = g_all[:p]  # (g_1, ..., g_p)
        g_last = g_all[p:]   # (g_{p+1}, ..., g_{p+q})

        # f(g_1, ..., g_p)
        f_idx_base = _multi_index(g_first, m) * r_m if p > 0 else 0
        f_values = f[f_idx_base:f_idx_base + r_m]

        # g(g_{p+1}, ..., g_{p+q})
        g_idx_base = _multi_index(g_last, m) * r_n if q > 0 else 0
        g_values = g[g_idx_base:g_idx_base + r_n]

        # Compute g_1 * g_2 * ... * g_p (product in G)
        prefix_product = group.identity_idx
        for gi in g_first:
            prefix_product = group.mult(prefix_product, gi)

        # Apply action of prefix_product to g_values
        acted_g = module_n.action(prefix_product) @ g_values

        # Tensor product: (f tensor acted_g)_{i,j} = f_i * acted_g_j
        tensor_values = np.outer(f_values, acted_g).flatten()

        # Place in result
        target_base = _multi_index(g_all, m) * r_tensor
        result[target_base:target_base + r_tensor] = tensor_values

    return result


def cup_product_matrix(
    group: FiniteGroup,
    module_m: GModule,
    module_n: GModule,
    p: int,
    q: int,
) -> np.ndarray:
    """
    Compute the cup product as a bilinear map (matrix form).

    Returns a matrix Cup of shape (dim C^{p+q}(G, M tensor N), dim C^p(G, M) * dim C^q(G, N))
    such that for f in C^p, g in C^q flattened:
        (f cup g) = Cup @ (f kron g)

    Actually, since cup product is bilinear, we return a 3-index tensor:
    Cup[out, f_idx, g_idx] gives the coefficient.
    For rank-1 trivial modules, this simplifies to a standard matrix.
    """
    m_ord = group.order
    r_m = module_m.rank
    r_n = module_n.rank
    r_tensor = r_m * r_n

    dim_f = r_m * (m_ord ** p)
    dim_g = r_n * (m_ord ** q)
    dim_out = r_tensor * (m_ord ** (p + q))

    # For rank-1 trivial modules, cup is a bilinear form
    # We compute it by evaluating on all basis pairs
    cup_tensor = np.zeros((dim_out, dim_f, dim_g), dtype=np.int64)

    for fi in range(dim_f):
        f_basis = np.zeros(dim_f, dtype=np.int64)
        f_basis[fi] = 1
        for gi in range(dim_g):
            g_basis = np.zeros(dim_g, dtype=np.int64)
            g_basis[gi] = 1
            result = cup_product_cochain(group, module_m, module_n, f_basis, g_basis, p, q)
            cup_tensor[:, fi, gi] = result

    return cup_tensor


def cup_product_trivial_z(
    group: FiniteGroup,
    f: np.ndarray,
    g: np.ndarray,
    p: int,
    q: int,
) -> np.ndarray:
    """
    Simplified cup product for trivial Z-coefficients (rank 1, trivial action).

    f in C^p(G, Z), g in C^q(G, Z) -> f cup g in C^{p+q}(G, Z).
    (f cup g)(g_1,...,g_{p+q}) = f(g_1,...,g_p) * g(g_{p+1},...,g_{p+q})
    """
    m = group.order

    result_dim = m ** (p + q)
    result = np.zeros(result_dim, dtype=np.int64)

    for target_tuple in cartesian_product(range(m), repeat=p + q):
        g_first = target_tuple[:p]
        g_last = target_tuple[p:]

        f_idx = _multi_index(g_first, m) if p > 0 else 0
        g_idx = _multi_index(g_last, m) if q > 0 else 0

        target_idx = _multi_index(target_tuple, m)
        result[target_idx] = int(f[f_idx]) * int(g[g_idx])

    return result


def verify_cup_leibniz(
    group: FiniteGroup,
    module_m: GModule,
    module_n: GModule,
    f: np.ndarray,
    g: np.ndarray,
    p: int,
    q: int,
) -> bool:
    """
    Verify the Leibniz rule: d(f cup g) = df cup g + (-1)^p f cup dg.

    This is the key property that makes cup product descend to cohomology.
    """
    from .cochains import coboundary_matrix

    # Build tensor product module
    tensor_module = _tensor_product_module(group, module_m, module_n)

    # d^{p+q}(f cup g)
    fg = cup_product_cochain(group, module_m, module_n, f, g, p, q)
    d_pq = coboundary_matrix(group, tensor_module, p + q)
    lhs = d_pq @ fg

    # d^p f cup g
    d_p = coboundary_matrix(group, module_m, p)
    df = d_p @ f
    df_cup_g = cup_product_cochain(group, module_m, module_n, df, g, p + 1, q)

    # (-1)^p * (f cup d^q g)
    d_q = coboundary_matrix(group, module_n, q)
    dg = d_q @ g
    f_cup_dg = cup_product_cochain(group, module_m, module_n, f, dg, p, q + 1)
    sign = (-1) ** p

    rhs = df_cup_g + sign * f_cup_dg

    return np.array_equal(lhs, rhs)


def _tensor_product_module(
    group: FiniteGroup, module_m: GModule, module_n: GModule
) -> GModule:
    """Construct M tensor N as a G-module with diagonal action."""
    r_m = module_m.rank
    r_n = module_n.rank
    r_tensor = r_m * r_n

    action_matrices = []
    for g in range(group.order):
        # Kronecker product of action matrices gives the tensor action
        tensor_action = np.kron(module_m.action(g), module_n.action(g))
        action_matrices.append(tensor_action)

    return GModule(group=group, rank=r_tensor, action_matrices=action_matrices)
