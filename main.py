import matplotlib.pyplot as plt
import numpy as np
import math


def func_d(x):
    return x / (x ** 2 + 1)


def func_i(x):
    return 1 / (x * np.sqrt(1 - np.log(x)))


# Левое разностное отношение
def left_diff_ratio(m, x, y):
    i = 1
    for j in range(1, len(x)):
        if x[j - 1] <= m < x[j]:
            i = j
            break
    return (y[i] - y[i - 1]) / (x[i] - x[i - 1])


# Правое разностное отношение
def right_diff_ratio(m, x, y):
    i = 0
    for j in range(len(x) - 1):
        if x[j] <= m < x[j + 1]:
            i = j
            break
    return (y[i + 1] - y[i]) / (x[i + 1] - x[i])


# Центральное разностное отношение
def central_diff_ratio(m, x, y):
    i = 0
    for j in range(1, len(x) - 1):
        if x[j - 1] <= m < x[j + 1]:
            i = j
            break
    return 0.5 * (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i])


#  Приближенное соотношение для второй производной
def second_der_approx_relation(m, x, y):
    i = 0
    for j in range(1, len(x) - 1):
        if x[j - 1] <= m < x[j + 1]:
            i = j
            break
    return (y[i - 1] - 2 * y[i] + y[i + 1]) / (x[i + 1] - x[i]) ** 2


# 1 производная, вычисленная по таблице производных
def der_func_d(x):
    return (-x ** 2 + 1) / (x ** 2 + 1) ** 2


# 2 производная, вычисленная по таблице производных
def der2_func_d(x):
    return (2 * x ** 3 - 6 * x) / (x ** 2 + 1) ** 3


# Правое разностное отношение для аналитических выражений
def an_right_diff_r(x, h):
    return (func_d(x + h) - func_d(x)) / h


# Левое разностное отношение для аналитических выражений
def an_left_diff_r(x, h):
    return (func_d(x) - func_d(x - h)) / h


# Центральное разностное отношение для аналитических выражений
def an_central_diff_r(x, h):
    return 0.5 * (func_d(x + h) - func_d(x - h)) / h


#  Приближенное соотношение для второй производной для аналитических выражений
def an_second_der_approx_relation(x, h):
    return (func_d(x - h) - 2 * func_d(x) + func_d(x + h)) / h ** 2


# Первообразная функции
def integral_f(x):
    return -2 * math.sqrt(1 - np.log(x))


#  Метод средних прямоугольников
def middle_rectangles_method(a, b, n):
    h = (b - a) / n
    integral = 0
    for i in range(1, n + 1):
        xi = a + h * (i - 0.5)
        integral += func_i(xi)
    return h * integral


def graph_middle_rectangles_method(a, b, h):
    y = []
    x = []
    xi = a
    while xi <= b:
        x.append(xi)
        y.append(func_i(xi))
        xi += h
    return x, y


#  Метод трапеций
def trapezium_method(a, b, n):
    integral = (func_i(a) + func_i(b)) / 2
    x = a
    h = (b - a) / n
    for i in range(1, n):
        x += h
        integral += func_i(x)
    return h * integral


#  Метод парабол
def parabola_method(a, b, n):
    if n % 2:
        n += 1
    integral = func_i(a) + func_i(b)
    x = a
    h = (b - a) / n
    for i in range(1, n):
        x += h
        if i % 2:
            integral += 4 * func_i(x)
        else:
            integral += 2 * func_i(x)
    return h / 3 * integral


def line_2_dotes(x1, x2, y1, y2, x):
    return (x - x1) * (y2 - y1) / (x2 - x1) + y1


def expand_list(in_arr, in_arg):
    arr = []
    arg = []
    for i in range(len(in_arr) - 1):
        x = np.arange(in_arg[i], in_arg[i + 1], (in_arg[i + 1] - in_arg[i]) / 100)
        for j in x:
            arg.append(j)
            arr.append(line_2_dotes(in_arg[i], in_arg[i + 1], in_arr[i], in_arr[i + 1], j))
    return arr, arg


def graph_dif1(x, h):
    plt.grid()
    plt.title("Graphs for 1 derivative, h = " + str(h))
    plt.plot(x, der_func_d(x), c='r', label='1 derivative calculated from the table of derivatives')
    plt.plot(x, an_right_diff_r(x, h), c='g', label='right difference relation')
    plt.plot(x, an_left_diff_r(x, h), c='m', label='left difference relation')
    plt.plot(x, an_central_diff_r(x, h), c='b', label='central difference relation')
    plt.legend()
    plt.show()
    plt.title("Error graph")
    plt.plot(x, abs(der_func_d(x) - an_central_diff_r(x, h)))
    plt.plot(x, abs(der_func_d(x) - an_left_diff_r(x, h)))
    plt.plot(x, abs(der_func_d(x) - an_right_diff_r(x, h)))
    plt.grid()
    plt.show()


def graph_dif2(x, h):
    plt.grid()
    plt.title("Graphs for 2 derivative, h = " + str(h))
    plt.plot(x, der2_func_d(x), c='r', label='2 derivative calculated from the table of derivatives')
    plt.plot(x, an_second_der_approx_relation(x, h), c='g', label='approximate relation for the second derivative')
    plt.legend()
    plt.show()
    plt.title("Error graph")
    plt.plot(x, abs(der2_func_d(x) - an_second_der_approx_relation(x, h)))
    plt.grid()
    plt.show()


def calc(power):
    n = 2
    while (abs(middle_rectangles_method(begin, end, n - 1) - middle_rectangles_method(begin, end, n))) > 10 ** power:
        n += 1
    h = (end - begin) / n
    print('\n10 to the ' + str(power) + ' power: ' + str(n))
    print("Result: " + str(round(middle_rectangles_method(begin, end, n), abs(power) + 1)))
    x, y = graph_middle_rectangles_method(begin, end, h)
    err, arg = expand_list(y, x)
    X = np.arange(begin, end, h / 100)

    plt.scatter(x, y, label="dote")
    plt.plot(X, func_i(X), label="original integral graph")
    plt.plot(arg, err, label='broken line error graph')
    plt.legend()
    plt.grid()
    plt.show()

    for i in range(len(arg)):
        err[i] = abs(func_i(arg[i]) - err[i])
    plt.plot(arg, err, label='error graph')
    plt.hlines(10 ** power, arg[0], arg[len(arg) - 1], color='r', label="graph 10 ^" + str(power), zorder=100)
    plt.legend()
    plt.grid()
    plt.show()


# 1.1.1
x = [0, 2, 4, 6]
y = [-2.5, 0.17, 2.5, 0.3]
a, b, c = 0.5, 2.5, 5.5
print('Left D.R.:')
print('value at point a: ' + str(left_diff_ratio(a, x, y)))
print('value at point b: ' + str(left_diff_ratio(b, x, y)))
print('value at point c: ' + str(left_diff_ratio(c, x, y)))
print('\nRight D.R.:')
print('value at point a: ' + str(right_diff_ratio(a, x, y)))
print('value at point b: ' + str(right_diff_ratio(b, x, y)))
print('value at point c: ' + str(right_diff_ratio(c, x, y)))
print('\nCenter:')
print('value at point a: ' + str(central_diff_ratio(a, x, y)))
print('value at point b: ' + str(central_diff_ratio(b, x, y)))
print('c: ' + str(central_diff_ratio(c, x, y)))
print('\nDiff. second level:')
print('value at point a: ' + str(second_der_approx_relation(a, x, y)))
print('value at point b: ' + str(second_der_approx_relation(b, x, y)))
print('value at point c: ' + str(second_der_approx_relation(c, x, y)))

#  1.1.2
h = [0.5, 0.1, 0.01]
X = np.arange(0, 10, 0.01)

# for diff 1
graph_dif1(X, h[0])
graph_dif1(X, h[1])
graph_dif1(X, h[2])

# for diff 2
graph_dif2(X, h[0])
graph_dif2(X, h[1])
graph_dif2(X, h[2])

# 1.2.1
begin, end = 1, np.exp(1) - 0.001
n = [1, 2, 4, 10, 50, 100]

print("\nIntegral value calculated by Newton's formula:")
print(integral_f(end) - integral_f(begin))

print("\nMedium rectangles method:")
for i in range(len(n)):
    print("n = " + str(n[i]))
    print(middle_rectangles_method(begin, end, n[i]))

print("\nTrapezium method:")
for i in range(len(n)):
    print("n = " + str(n[i]))
    print(trapezium_method(begin, end, n[i]))

print("\nParabola method:")
for i in range(len(n)):
    print("n = " + str(n[i]))
    print(parabola_method(begin, end, n[i]))

# 1.2.2
calc(-1)
calc(-3)
calc(-5)
calc(-8)

# 1.2.3
n = 100
a = 0.5
b = 2.5
s1 = trapezium_method(a, b, n)
s2 = integral_f(b) - integral_f(a)
print("\nExact value: " + str(s2))
print('Method of trapezoid: ' + str(s1))
