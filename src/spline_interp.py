import src.equation_systems.equation_systems as eqs

import numpy as np


def get_equation(X_red, Y_red):
    # Creating the matrices M and b such that
    # "M x spl_coeefs = B"
    # is the equation that can be solved for spl_coeffs.

    points_n = X_red.shape[0]
    splines_n = points_n - 1
    M = []
    B = []

    # (I) S_i(x_i) = y_i
    for i in range(splines_n):
        # 1 for a_i, 0 for rest of the splines coefficients
        M.append([
            1 if j == 4 * i else 0
            for j in range(4 * splines_n)
        ])
        B.append([Y_red[i]])

    # (II) S_i(x_{i+1}) = y_{i+1}
    for i in range(splines_n):
        # 1 for a_i, x_{i+1}-x_i for b_i, (x_{i+1}-x_i)^2 for c_i, (x_{i+1}-x_i)^3 for d_i
        M.append([
            1 if j == 4 * i else
            X_red[i + 1] - X_red[i] if j == 4 * i + 1 else
            (X_red[i + 1] - X_red[i])**2 if j == 4 * i + 2 else
            (X_red[i + 1] - X_red[i])**3 if j == 4 * i + 3 else 0
            for j in range(4 * splines_n)
        ])
        B.append([Y_red[i + 1]])

    # (III) S_j'(x_{j+1}) = S_{j+1}'(x_{j+1})
    for i in range(splines_n - 1):
        # 1 for b_i, -1 for b_{i+1}, 2*(x_{i+1}-x_i) for c_i, 3*(x_{i+1}-x_i)^2 for d_i
        M.append([
            1 if j == 4 * i + 1 else
            2 * (X_red[i + 1] - X_red[i]) if j == 4 * i + 2 else
            3 * (X_red[i + 1] - X_red[i])**2 if j == 4 * i + 3 else
            -1 if j == 4 * (i + 1) + 1 else 0
            for j in range(4 * splines_n)
        ])
        B.append([0])

    # (IV) S_j''(x_{j+1}) = S_{j+1}''(x_{j+1})
    for i in range(splines_n - 1):
        # 2 for c_i, 6*(x_{i+1}-x_i) for d_i, -2 for c_{i+1}
        M.append([
            2 if j == 4 * i + 2 else
            6 * (X_red[i + 1] - X_red[i]) if j == 4 * i + 3 else
            -2 if j == 4 * (i + 1) + 2 else 0
            for j in range(4 * splines_n)
        ])
        B.append([0])

    # (V) S_0''(x_0) = 0
    # 1 for c_0
    M.append([
        1 if j == 2 else 0
        for j in range(4 * splines_n)
    ])
    B.append([0])
    # S_{n-2}''(x_{n-1}) = 0
    # 1 for c_{n-2}, 3*(x_{n-1}-x_{n-2}) for d_{n-2}
    M.append([
        1 if j == (splines_n - 1) * 4 + 2 else
        3 * (X_red[splines_n - 1] - X_red[splines_n - 2]) if j ==  (splines_n - 1) * 4 + 3 else 0
        for j in range(4 * splines_n)
    ])
    B.append([0])

    return M, B


def get_interp_values(X_interp, X_red, Y_red):
    M, B = get_equation(X_red, Y_red)
    spl_coeffs = np.linalg.solve(M, B)
    Y_interp = []

    for i in range(X_red.shape[0] - 1):
        for x in range(round(X_red[i]), round(X_red[i + 1] +
                                              (1 if i == X_red.shape[0] - 2 else 0))):
            Y_interp.append(
                spl_coeffs[4 * i][0] +
                spl_coeffs[4 * i + 1][0] * (x - X_red[i]) +
                spl_coeffs[4 * i + 2][0] * (x - X_red[i])**2 +
                spl_coeffs[4 * i + 3][0] *  (x - X_red[i])**3
            )

    return Y_interp
