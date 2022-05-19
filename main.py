import sympy as sp


def first_task():
    print('\nFirst Task\n')
    x = [0.43, 0.48, 0.55, 0.62, 0.70, 0.75]
    y = [1.63597, 1.73234, 1.87686, 2.03345, 2.22846, 2.35973]
    findx = [0.702, 0.512, 0.645, 0.736, 0.608]

    # Lagrange polynomial

    xsym = sp.Symbol('x')
    polynomial = 0  # F
    for i in range(len(x)):
        y_mul = y[i]
        fraction = 1
        for j in range(len(x)):
            if i != j:
                fraction *= (xsym - x[j]) / (x[i] - x[j])
        polynomial += y_mul * fraction

    polynomial = sp.expand(polynomial)
    print(polynomial)
    sp.plot(polynomial, (xsym, 0, 1))

    for i in range(len(findx)):
        print(polynomial.subs(xsym, findx[i]))


def newton_interpolation(x: float):
    table = [[1.415,
              1.420,
              1.425,
              1.430,
              1.435,
              1.440,
              1.445,
              1.450,
              1.455,
              1.460,
              1.465],
             [0.888551,
              0.889599,
              0.890637,
              0.891667,
              0.892687,
              0.893698,
              0.894700,
              0.895693,
              0.896677,
              0.897653,
              0.898619]]

    # All deltas finder:

    for i in range(0, len(table[1]) - 1):
        table.append([])

    for i in range(2, len(table)):
        for j in range(1, len(table[i - 1])):
            table[i].append(table[i - 1][j] - table[i - 1][j - 1])

    # Find x0:
    # Default case: y0 is table[2][0], where 0 is default index
    # min_sub is x - x0, but x0 is to find, so it is x - table[0][0] by default
    default_index = 0
    min_sub = abs(x - table[0][default_index])
    x_zero = table[0][default_index]

    # In the next loop we are looking for the nearest value for x - x0
    for i in range(1, len(table[0])):
        if abs(x - table[0][i]) < min_sub:
            min_sub = x - table[0][i]
            x_zero = table[0][i]
            default_index = i
    # print(x_zero)

    # t = (x-x0)/h
    h = abs(table[0][1] - table[0][0])
    t = (x - x_zero) / h

    # y0
    y_zero = table[1][default_index]

    # Now we have all data for calculating newton interpolation
    # Forward and backward interpolation cases:

    y = 0
    pre_y_modifier = 1
    # print(f'{x} |{x_zero} | {y_zero} | {h} | {t}')
    if x < table[0][(len(table[0]) - 1) // 2]:
        for i in range(0, (len(table[default_index + 2]))):
            y += y_zero * pre_y_modifier
            pre_y_modifier *= ((t - i) * table[i + 2][default_index]) / (i + 1)
    else:
        y = y_zero
        for i in range(0, len(table) - 2):
            pre_y_modifier *= (t - i)*table[i+2][len(table[i + 2]) - 1]/(i + 1)
            y += y_zero * pre_y_modifier

    return y


def second_task():
    values = [1.4161, 1.4625, 1.4135, 1.470]
    print('\nSecond Task\n')
    for value in values:
        print(f'y({value}) = {newton_interpolation(value)}')


if __name__ == '__main__':
    first_task()
    second_task()
