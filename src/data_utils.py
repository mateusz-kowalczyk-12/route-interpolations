import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# def get_reduced_data(X, Y, points_n):
#     center = int(X.shape[0] / 2)
#     step = (2 * (X.shape[0] - 1)) / ((1 + int(points_n / 2)) * points_n)
#
#     red_idxs = [center + i * step
#                 for i in range(int(points_n / 2 + 1))]
#     lower_idxs = [center - i * step
#                   for i in range(1, int(points_n / 2))]
#     lower_idxs.reverse()
#
#     red_idxs = lower_idxs + red_idxs
#
#     X_red, Y_red = X[red_idxs], Y[red_idxs]
#     X_red = np.concatenate([X_red, [X[X.shape[0] - 1]]])
#     Y_red = np.concatenate([Y_red, [Y[Y.shape[0] - 1]]])
#
#     return X_red, Y_red


def get_reduced_data(X, Y, points_n):
    step = X.shape[0] / (points_n - 1)
    red_idxs = [int(i * step)
                for i in range(points_n - 1)]

    X_red, Y_red = X[red_idxs], Y[red_idxs]
    X_red = np.concatenate([X_red, [X[X.shape[0] - 1]]])
    Y_red = np.concatenate([Y_red, [Y[Y.shape[0] - 1]]])

    return X_red, Y_red


def get_data(filename, sep, points_n):
    dataset = pd.read_csv(f'datasets/{filename}', sep=sep, header=None, names=['x', 'y'])
    X, Y = dataset['x'].values, dataset['y'].values
    X_red, Y_red = get_reduced_data(X, Y, points_n)

    return X, Y, X_red, Y_red


def plot(X, Y, X_red, Y_red, X_interp, Y_interp, ylim, method, points_n, route):
    method_name = 'Lagrange\'a' if method == 'L' else 'funkcjami sklejanymi'
    title = f'Interpolacja {method_name} profilu wysokościowego trasy "{route}"'

    plt.figure(figsize=(16, 10))
    plt.plot(X_interp, Y_interp, label='funkcja interpolująca', zorder=1)
    plt.plot(X, Y, label='funkcja bazowa', zorder=0)
    plt.scatter(X_red, Y_red, s=20, c='darkred', zorder=2, label=f'węzły interpolacji ({points_n})')

    plt.ylim(ylim)

    plt.title(title, fontsize=18)
    plt.xlabel('odległość [m]', fontsize=16)
    plt.ylabel('wysokość [m n.p.m.]', fontsize=16)

    plt.legend(fontsize=14)

    plt.savefig(f'report/media/{route}_{method}_{points_n}.png', bbox_inches='tight')
    plt.show()
