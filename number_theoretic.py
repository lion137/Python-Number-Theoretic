# Chinese Remainder Theorem algorithm
# https://en.wikipedia.org/wiki/Chinese_remainder_theorem

# Modular exponeneation also here: https://github.com/lion137/Python-Algorithms/blob/master/modular_exponenation.py

# Legendre / Jacoby Symbols algorithms
# https://en.wikipedia.org/wiki/Jacobi_symbol
# https://en.wikipedia.org/wiki/Legendre_symbol

# Quadratic Residue Algorithm
# https://en.wikipedia.org/wiki/Quadratic_residue


from functools import reduce
from operator import add
from Euler import modularexponenation, exp
from random import randint
from math import log2


def chinese_remainder(p, n):
    """takes vector p of coprimes and vector n of residues,
    returns solution to Chinese Remainder Theorem relations"""

    # helper function
    def mul_inv(a, n):
        """returns modular multiplicative inverse of a mod n"""
        for x in range(1, 10 * 100):
            if x * a % n == 1:
                return x

    m = reduce(lambda x, y: x * y, p)
    mi = list(map(lambda x: m // x, p))
    v_inverses = [mul_inv(x, y) for [x, y] in zip(mi, p)]
    return reduce(add, [x * y * z for [x, y, z] in zip(n, v_inverses, mi)]) % m


def legendre_jacoby(m, a):
    """takes odd integer m and integer a,
    returns Jacoby Symbol (a/m), which if
    m is prime not = 2 is Legendre Symbol too"""
    even = lambda x: x % 2 == 0
    a = a % m
    t = 1
    while a != 0:
        while even(a):
            a //= 2
            if 3 <= m % 8 <= 5:
                t = -t
        a, m = m, a
        if a == m == 3 % 4:
            t = -t
        a = a % m
    if m == 1:
        return t
    return 0


def square_roots_mod_p(p, a):
    """takes prime p > 2 and integer a such that (a/p) = 1 (Legendre Symbol)
    returns solution x to x^2 = a (mod p) (quadratic residue)"""

    # Simple p = 3, 5, 7 mod 8:
    a = a % p
    if p % 8 == 3 or p % 8 == 7:
        x = modularexponenation(a, (p + 1) // 4, p)
        return x
    if p % 8 == 5:
        x = modularexponenation(a, (p + 3) // 8, p)
        c = x ** 2 % p
        if c % a != p:
            x = (x * modularexponenation(2, (p - 1) // 4, p)) % p
        return x
    if p % 8 == 1:
        for _ in range(10 ** 100):
            d = randint(2, p - 1)
            if legendre_jacoby(p, d) == -1:
                break
        loop_cnt = 0
        t = 0
        for s0 in range(int(log2(p)), 0, -1):
            tmp = exp(2, s0)
            for t0 in range(1, (p - 1) // 2 + 1, 2):
                if p - 1 == tmp * t0:
                    loop_cnt = 1
                    t = t0
                    s = s0
                    break
            if loop_cnt == 1:
                break
        A = (a ** t) % p
        D = (d ** t) % p
        m = 0
        for x in range(10**10):
            if (A * D**x) % p == 1:
                m = x
                break
        x = (a ** ((t + 1) // 2)) * D**(m // 2) % p
        return x


if __name__ == '__main__':
    p = [3, 5, 7]
    a = [2, 2, 6]
    a1 = [2, 3, 2]
    #q = 65
    # print(chinese_remainder(p, a)) # -> 62 Example from the Number Theory textbook.
    print(legendre_jacoby(11, 15))
    print(square_roots_mod_p(11, 15))
    print(81 % 11)
    '''        for i in range(s):
            if ((A * D**m)**(2**(s - 1 - i))) % p == -1:
                m += exp(2, i)'''




