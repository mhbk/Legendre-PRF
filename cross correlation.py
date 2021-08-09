from sympy.ntheory import legendre_symbol
import matplotlib.pyplot as plt
import math
import numpy as np


def jacobi(a, n):
    if a == 0:
        return 0
    assert (n > a > 0 and n % 2 == 1)
    t = 1
    while a != 0:
        while a % 2 == 0:
            a /= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        a %= n
    if n == 1:
        return t
    else:
        return 0


def compute_legendre_arr(p: int):
    arr = np.zeros(p, dtype=int)
    for i in range(0, p):
        arr[i] = jacobi(i, p)
    return arr


prime = []
N = 1000
cross_correlation = 0
L = legendre_symbol
y = []
x = []

y_r = []

for i in range(2, N):
    sq_i = int(math.sqrt(i)) + 1
    temp = True
    for j in prime:
        if j > sq_i:
            break
        if i % j == 0:
            temp = False
            break
    if temp:
        prime.append(i)

for p in prime[2:-1]:
    print(p)
    arr_random = np.random.choice([0, 1], size=(p,), p=[1. / 2, 1. / 2])
    arr_Legendre = compute_legendre_arr(p)
    cross_correlation_ave = 0
    cross_correlation_ave_rand = 0
    for i in range(1, p):
        cross_correlation = 0
        cross_correlation_rand = 0
        for j in range(0, p - 1):
            if (arr_Legendre[(j + i) % p] - arr_Legendre[j % p]) % 2 == 0:
                cross_correlation += 1
            else:
                cross_correlation += -1
            # cross_correlation += (-1) ** (jacobi(j + i, prime) - jacobi(j, prime))
            if (arr_random[(j + i) % p] - arr_random[j % p]) % 2 == 0:
                cross_correlation_rand += 1
            else:
                cross_correlation_rand += -1
        cross_correlation_ave += cross_correlation
        cross_correlation_ave_rand += cross_correlation_rand
    cross_correlation_ave = cross_correlation_ave / (p - 1)
    cross_correlation_ave_rand = cross_correlation_ave_rand / (p - 1)
    y.append(cross_correlation_ave)
    y_r.append(cross_correlation_ave_rand)
    x.append(p)

plt.plot(x, y, label="Legendre PRF")
plt.plot(x, y_r, label="random sequence")
plt.legend()
plt.show()
