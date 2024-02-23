import copy
import math


def build_Ab(a1, a2, a3, f, N):
    A = []
    b = []

    for i in range(N):
        b.append([math.sin(i * (f + 1))])
        A.append([])

        for j in range(N):
            element_to_append = None
            if j == i:
                element_to_append = a1
            elif abs(i - j) == 1:
                element_to_append = a2
            elif abs(i - j) == 2:
                element_to_append = a3
            else:
                element_to_append = 0

            A[i].append(element_to_append)

    return A, b


def build_ULD(A):
    U = copy.deepcopy(A)
    L = copy.deepcopy(A)
    D = copy.deepcopy(A)

    for row in range(len(A)):  # for each row
        for col in range(len(A[0])):  # for each column
            if row < col:
                U[row][col] = -U[row][col]
                L[row][col] = 0
                D[row][col] = 0
            elif row == col:
                U[row][col] = 0
                L[row][col] = 0
            elif row > col:
                U[row][col] = 0
                L[row][col] = -L[row][col]
                D[row][col] = 0

    return U, L, D


def LU_factorise(A):
    """Returns the matrices L and U, such that A = LU"""

    L = build_I(len(A))
    U = copy.deepcopy(A)

    for k in range(len(A) - 1):
        for j in range(k + 1, len(A)):
            L[j][k] = U[j][k] / U[k][k]

            for i in range(k, len(A)):
                U[j][i] -= L[j][k] * U[k][i]

    return L, U


def build_I(m):
    """Returns an identity matrix I of dimensions m x m"""

    I = []

    for i in range(m):
        I.append([])

        for j in range(m):
            I[i].append(1 if i == j else 0)

    return I


def add_matrices(A, B):
    """Returns the sum C = A + B"""

    if (not len(A) == len(B)) or (not len(A[0]) == len(B[0])):
        raise Exception('Trying to add matrices with not equal dimensions')

    C = copy.deepcopy(A)

    for row in range(len(A)):
        for col in range(len(A[0])):
            C[row][col] += B[row][col]

    return C


def subtract_matrices(A, B):
    """Returns the difference C = A - B"""

    minus_B = copy.deepcopy(B)

    for i in range(len(B)):
        for j in range(len(B[0])):
            minus_B[i][j] = -minus_B[i][j]

    return add_matrices(A, minus_B)


def forward_substitution(A, b):
    """A is a lower triangular matrix, b is a vector.
    Returns the solution x of the equation Ax = b"""

    x = [[b[0][0] / A[0][0]]]

    for i in range(1, len(A)):  # for each row
        row_sum = 0
        for j in range(i):
            row_sum += x[j][0] * A[i][j]

        x.append([(b[i][0] - row_sum) / A[i][i]])

    return x


def backward_substitution(A, b):
    """A is an upper triangular matrix, b is a vector.
    Returns the solution x of the equation Ax = b"""

    x = [[b[len(A) - 1][0] / A[len(A) - 1][len(A) - 1]]]

    for i in range(len(A) - 2, -1, -1):  # for each row
        row_sum = 0
        for j in range(len(A) - 1, i, -1):
            row_sum += x[j - i - 1][0] * A[i][j]

        x.insert(0, [(b[i][0] - row_sum) / A[i][i]])

    return x


def diag_inv(D):
    """Returns the inversion of the diagonal matrix D"""

    D_inv = copy.deepcopy(D)

    for i in range(len(D_inv)):
        D_inv[i][i] = 1 / D_inv[i][i]

    return D_inv


def multiply_matrices(A, B):
    """Returns the product C = AB"""

    if not len(A[0]) == len(B):
        raise Exception('Trying to multiply not compatible matrices')

    p = len(A)  # numbers of rows of A
    q = len(B)  # numbers of columns of A/number of rows of B
    r = len(B[0])  # number of columns of B

    C = []  # will be of dimensions p x r

    for i in range(p):  # for each row of the result matrix
        C.append([])

        for j in range(r):  # for each column of the result matrix
            new_val = 0
            for k in range(q):
                new_val += A[i][k] * B[k][j]

            C[i].append(new_val)

    return C


def multiply_diag_by_matrix(D, A):
    """D is a diagonal matrix
    Returns the product C = DA"""

    if not len(D[0]) == len(A):
        raise Exception('Trying to multiply not compatible matrices')

    p = len(D)  # numbers of rows of D
    q = len(A)  # numbers of columns of A/number of rows of A
    r = len(A[0])  # number of columns of A

    C = []  # will be of dimensions p x r

    for i in range(p):  # for each row of the result matrix
        C.append([])

        for j in range(r):  # for each column of the result matrix
            C[i].append(D[i][i] * A[i][j])

    return C


def residual_error(A, x, b):
    """Returns the difference Ax - b"""

    return subtract_matrices(multiply_matrices(A, x), b)


def norm2(v):
    """Returns the 2. norm of the vector v"""

    norm = 0
    for vi in v:
        norm += vi[0] ** 2

    return math.sqrt(norm)
