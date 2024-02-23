import src.equation_systems.data_tools as dt
import math
import numpy as np
import time


def jacobi_solve(A, b):
    prep_time = time.time()

    U, L, D = dt.build_ULD(A)
    D_inv_L_plus_U = dt.multiply_diag_by_matrix(dt.diag_inv(D), dt.add_matrices(L, U))
    D_inv_b = dt.forward_substitution(D, b)

    x = [[1] for _ in range(len(b))]
    res = dt.residual_error(A, x, b)
    res_norm = dt.norm2(res)

    prep_time = time.time() - prep_time

    iteration = 0
    res_norms = [res_norm]

    iter_time = time.time()

    while not math.isinf(res_norm) and res_norm > 10**-9 and iteration < 600:
        iteration += 1

        x = dt.add_matrices(dt.multiply_matrices(D_inv_L_plus_U, x), D_inv_b)
        res = dt.residual_error(A, x, b)
        res_norm = dt.norm2(res)

        res_norms.append(res_norm)

    iter_time = time.time() - iter_time

    return x, iteration, res_norms, prep_time, iter_time


def gauss_seidl_solve(A, b):
    prep_time = time.time()

    U, L, D = dt.build_ULD(A)
    D_minus_L = dt.subtract_matrices(D, L)
    D_minus_L_inv_b = dt.forward_substitution(D_minus_L, b)

    x = [[1] for _ in range(len(b))]
    res = dt.residual_error(A, x, b)
    res_norm = dt.norm2(res)

    prep_time = time.time() - prep_time

    iteration = 0
    res_norms = [res_norm]

    iter_time = time.time()

    while not math.isinf(res_norm) and res_norm > 10 ** -9 and iteration < 600:
        iteration += 1

        x = dt.add_matrices(dt.forward_substitution(D_minus_L, dt.multiply_matrices(U, x)), D_minus_L_inv_b)
        res = dt.residual_error(A, x, b)
        res_norm = dt.norm2(res)

        res_norms.append(res_norm)

    iter_time = time.time() - iter_time

    return x, iteration, res_norms, prep_time, iter_time


def LU_fact_solve(A, b):
    fact_time = time.time()
    L, U = dt.LU_factorise(A)
    fact_time = time.time() - fact_time

    subst_time = time.time()
    y = dt.forward_substitution(L, b)
    x = dt.backward_substitution(U, y)
    subst_time = time.time() - subst_time

    res = dt.residual_error(A, x, b)
    res_norm = dt.norm2(res)

    return x, res_norm, fact_time, subst_time

