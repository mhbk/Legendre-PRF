EPSILON = 0.0001


# Recursively calculates the gcd of two polynomials in given finite field p (for prime p)
# Polynomials are given by a list of coefficients from largest to smallest.
# When p=0 tries to calculate the gcd in R, percision makes this difficult, and is not reliable.
def gcd(f, g, p=0, verbose=False):
    if len(f) < len(g):
        return gcd(g, f, p, verbose)

    r = [0] * len(f)
    r_mult = reciprocal(g[0], p) * f[0]

    for i in range(len(f)):
        if i < len(g):
            r[i] = f[i] - g[i] * r_mult
        else:
            r[i] = f[i]
        if (p != 0):
            r[i] %= p

    if verbose:
        print(f, 'by', g, 'got', r)

    while (abs(r[0]) < EPSILON):
        r.pop(0)
        if (len(r) == 0):
            return g

    return gcd(r, g, p, verbose)


# returns reciprocal of n in finite field of prime p, if p=0 returns 1/n#
def reciprocal(n, p=0):
    if p == 0:
        return 1 / n
    for i in range(p):
        if (n * i) % p == 1:
            return i
    return None

def gcd_modify(f, g, p=0, verbose=False):
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



f = [-1,1,0]
g = [1,0,-1]
print(gcd(f, g, 2))
print(gcd_modify(f, g, 2))
# print(reciprocal(12, 19))
