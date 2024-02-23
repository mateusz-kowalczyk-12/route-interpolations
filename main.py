import src.data_utils as du
import src.lagrange_interp as lagr
import src.spline_interp as spl


def plain_terrain(method, points_n):
    X, Y, X_red, Y_red = du.get_data('chelm.csv', ' ', points_n)
    X_interp = [x for x in range(round(X_red[X_red.shape[0] - 1]) + 1)]

    if method == 'L':
        Y_interp = lagr.get_interp_values(X_interp, X_red, Y_red)
    else:
        Y_interp = spl.get_interp_values(X_interp, X_red, Y_red)

    du.plot(X, Y, X_red, Y_red, X_interp, Y_interp, [0, 100], method, points_n, 'Chełm')


def single_depth(method, points_n):
    X, Y, X_red, Y_red = du.get_data('GlebiaChallengera.csv', ',', points_n)
    X_interp = [x for x in range(round(X_red[X_red.shape[0] - 1]) + 1)]

    if method == 'L':
        Y_interp = lagr.get_interp_values(X_interp, X_red, Y_red)
    else:
        Y_interp = spl.get_interp_values(X_interp, X_red, Y_red)

    du.plot(X, Y, X_red, Y_red, X_interp, Y_interp, [-12000, 0], method, points_n, 'Głębia Challengera')


def many_heights(method, points_n):
    X, Y, X_red, Y_red = du.get_data('rozne_wniesienia.csv', ',', points_n)
    X_interp = [x for x in range(round(X_red[X_red.shape[0] - 1]) + 1)]

    if method == 'L':
        Y_interp = lagr.get_interp_values(X_interp, X_red, Y_red)
    else:
        Y_interp = spl.get_interp_values(X_interp, X_red, Y_red)

    du.plot(X, Y, X_red, Y_red, X_interp, Y_interp, [0, 1200], method, points_n, 'Różne wzniesienia')


if __name__ == '__main__':
    while(True):
        terrain = input('1 - chelm, 2 - GlebiaChallengera, 3 - rozne_wniesienia\n')
        method = input('L - Lagrange method, S - splines method\n')
        points_n = int(input('points number: '))
        match(terrain):
            case '1':
                plain_terrain(method, points_n)
            case '2':
                single_depth(method, points_n)
            case '3':
                many_heights(method, points_n)
