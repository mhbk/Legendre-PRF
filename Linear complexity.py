import math
import sys

import matplotlib.pyplot as plt
import numpy
import copy
import scipy.special as spc


def linear_complexity(self, bin_data, block_size=500):
    """
    Note that this description is taken from the NIST documentation [1]
    [1] http://csrc.nist.gov/publications/nistpubs/800-22-rev1a/SP800-22rev1a.pdf
    The focus of this test is the length of a linear feedback shift register (LFSR). The purpose of this test is to
    determine whether or not the sequence is complex enough to be considered random. Random sequences are
    characterized by longer LFSRs. An LFSR that is too short implies non-randomness.
    :param bin_data: a binary string
    :param block_size: the size of the blocks to divide bin_data into. Recommended block_size >= 500
    :return:
    """
    dof = 6
    piks = [0.01047, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833]

    t2 = (block_size / 3.0 + 2.0 / 9) / 2 ** block_size
    mean = 0.5 * block_size + (1.0 / 36) * (9 + (-1) ** (block_size + 1)) - t2

    num_blocks = int(len(bin_data) / block_size)
    if num_blocks > 1:
        block_end = block_size
        block_start = 0
        blocks = []
        for i in range(num_blocks):
            blocks.append(bin_data[block_start:block_end])
            block_start += block_size
            block_end += block_size

        complexities = []
        for block in blocks:
            complexities.append(self.berlekamp_massey_algorithm(block))

        t = ([-1.0 * (((-1) ** block_size) * (chunk - mean) + 2.0 / 9) for chunk in complexities])
        vg = numpy.histogram(t, bins=[-9999999999, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 9999999999])[0][::-1]
        im = ([((vg[ii] - num_blocks * piks[ii]) ** 2) / (num_blocks * piks[ii]) for ii in range(7)])

        chi_squared = 0.0
        for i in range(len(piks)):
            chi_squared += im[i]
        p_val = spc.gammaincc(dof / 2.0, chi_squared / 2.0)
        return p_val
    else:
        return -1.0


def berlekamp_massey_algorithm(block_data):
    """
    An implementation of the Berlekamp Massey Algorithm. Taken from Wikipedia [1]
    [1] - https://en.wikipedia.org/wiki/Berlekamp-Massey_algorithm
    The Berlekamp–Massey algorithm is an algorithm that will find the shortest linear feedback shift register (LFSR)
    for a given binary output sequence. The algorithm will also find the minimal polynomial of a linearly recurrent
    sequence in an arbitrary field. The field requirement means that the Berlekamp–Massey algorithm requires all
    non-zero elements to have a multiplicative inverse.
    :param block_data:
    :return:
    """
    n = len(block_data)
    c = numpy.zeros(n)
    b = numpy.zeros(n)
    c[0], b[0] = 1, 1
    l, m, i = 0, -1, 0
    int_data = [int(el) for el in block_data]
    while i < n:
        v = int_data[(i - l):i]
        v = v[::-1]
        cc = c[1:l + 1]
        d = (int_data[i] + numpy.dot(v, cc)) % 2
        if d == 1:
            temp = copy.copy(c)
            p = numpy.zeros(n)
            for j in range(0, l):
                if b[j] == 1:
                    p[j + i - m] = 1
            c = (c + p) % 2
            if l <= 0.5 * i:
                l = i + 1 - l
                m = i
                b = temp
        i += 1
    return l


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
    arr = numpy.zeros(p, dtype=int)
    for i in range(0, p):
        arr[i] = jacobi(i, p)
    return arr


def compute_h_s_x(sequences):
    length = len(sequences)
    index = length
    for i in range(length - 1, 0, -1):
        if sequences[i] == 0:
            index -= 1
        else:
            break
    ans = numpy.zeros(index, dtype=int)
    for i in range(0, index):
        ans[i] = sequences[index - 1 - i]
    return ans

def compute_x_n(degree:int):
    arr = numpy.zeros(degree, dtype=int)
    arr[0] = 1
    arr[-1] = -1
    return arr


def prime_calculator(upper_bound):
    prime = []

    for i in range(2, upper_bound):
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
    return prime

EPSILON = 0.0001


# Recursively calculates the gcd of two polynomials in given finite field p (for prime p)
# Polynomials are given by a list of coefficients from largest to smallest.
# When p=0 tries to calculate the gcd in R, percision makes this difficult, and is not reliable.
def gcd(f, g, p=0, verbose=False):
    if len(f) < len(g):
        return gcd(g, f, p, verbose)
    while True:
        r = [0] * len(f)
        r_mult = reciprocal(g[0], p) * f[0]
        for i in range(len(f)):
            if i < len(g):
                r[i] = f[i] - g[i] * r_mult
            else:
                r[i] = f[i]
            if p != 0:
                r[i] %= p

        while abs(r[0]) < EPSILON:
            r.pop(0)
            if len(r) == 0:
                return g
        f = g
        g = r

    # return gcd(r, g, p, verbose)


# returns reciprocal of n in finite field of prime p, if p=0 returns 1/n#
def reciprocal(n, p=0):
    if n==0:
        return 0
    if p == 0:
        return 1 / n
    for i in range(0,p):
        if (n * i) % p == 1:
            return i
    return 0
# sys.setrecursionlimit(10000)
# arr = [0, 1, 0, 1, 1]
# print(berlekamp_massey_algorithm(arr))
N = 1000
x = []
y_Legendre = []
y_random = []
y_divided = []
prime = prime_calculator(N)
for p in prime[1:-1]:
    # bound = min(30 , int((p-1)/2))
    bound = int((p-1)/2)
    Legendre_arr = compute_legendre_arr(p)[0:bound]
    arr_random = numpy.random.choice([0, 1], size=(p,), p=[1. / 2, 1. / 2])[0:bound]
    poly1 = compute_x_n(bound)
    poly2 = compute_h_s_x(Legendre_arr)
    poly3 = compute_h_s_x(arr_random)
    temp1 = gcd(Legendre_arr, poly1,p)
    temp2 = gcd(poly3, poly1 , p)
    # print(temp1)
    # print("")
    linear_complex_legendre = bound - len(temp1) -1
    linear_complex_random = bound - len(temp2) - 1
    y_Legendre.append(linear_complex_legendre)
    y_random.append(linear_complex_random)
    y_divided.append(p/2)
    x.append(p)

plt.plot(x, y_Legendre, label="Legendre")
plt.plot(x, y_random, label="Random")
plt.plot(x, y_divided, label="y = x/2")
plt.legend()
plt.show()

