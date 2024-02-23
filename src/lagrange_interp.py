import numpy as np


def get_Fi(i, x, X_red):
    value = 1
    for j in range(X_red.shape[0]):
        if not i == j:
            value *= (x - X_red[j]) / (X_red[i] - X_red[j])
    return value


def get_interp_value(x, X_red, Y_red):
    value = 0
    for i in range(X_red.shape[0]):
        value += Y_red[i] * get_Fi(i, x, X_red)
    return value


def get_interp_values(X_interp, X_red, Y_red):
    Y_interp = []
    for x in X_interp:
        if x % 1000 == 0:
            print(f'{x}/{len(X_interp)}')
        Y_interp.append(get_interp_value(x, X_red, Y_red))
    return Y_interp
