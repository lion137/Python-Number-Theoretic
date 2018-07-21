# Factorization algorithms

# 1. Simple trial, having list of enough primes to factorize certain integer
# Factorize by trial division over the list, we assume, that we have, a list of
# primes generated, for example from Sieve of Erasthotenes.

# Function factorize_trial

from Euler import miller_rabin, prime_sieve, modularexponenation
import random
from functools import reduce
import operator as op
from math import gcd, sqrt, fabs, log2, log

def factorize_trial(n, primes):
    """Factorizing n"""
    if miller_rabin(n):
        return n
    result = []
    for p in primes:
        while n % p == 0:
            result.append(p)
            n //= p
        if miller_rabin(n):
            result.append(n)
            break
    return result


def PollardRho(n):
    if miller_rabin(n):
        return n
    out = 0
    i = 0
    xi = random.randint(0,n-1)
    xi = 5
    k = 2
    y = xi
    while True:
        i = i + 1
        xi = ((xi ** 2) + 1) % n
        d = gcd(y - xi,n)
        if d != 1 and d != n:
            if miller_rabin(d):
                out = d
                break
        if i == k:
            y = xi
            k = 2 * k
    return out

def rho(n):
    if miller_rabin(n):
        return n
    def g(x, m):
        return (x*x - 1) % m
    x = 2
    y = 2
    d = 1
    while d == 1:
        x = g(x, n)
        y = g(g(y, n), n)
        d = gcd(abs(x - y), n)
    if d == n:
        return False
    else:
        return d


def fermat_factorization(n):
    def is_square(a):
        if a < 0:
            a = -a
        return int(sqrt(a)) * int(sqrt(a)) == a
    a = int(sqrt(n))
    b = a * a - n
    while not is_square(b):
        a += 1
        b = a * a - n
    return [a - int(sqrt(b)), a + int(sqrt(b))]


def pollard_p_minus_1(n, max_iter, b):
    """Performing Pollard p- 1 factorization,
    with odd number input"""
    for _ in range(max_iter):
        primes = prime_sieve(b)
        for i in range(len(primes)):
            primes[i] = primes[i] ** (int(log(n, primes[i])))
        m = reduce(op.mul, primes)
        g = gcd(modularexponenation(2, m, n) - 1, n)
        if 1 < g < n:
            return g
        if g == 1:
            b += (b // 5)
        if g == n:
            b -= (b // 5)
'''function gcd(a, b)
    while a ≠ b 
        if a > b
           a := a − b; 
        else
           b := b − a; 
    return a;'''


if __name__ == '__main__':

    import time

    start = time.time()
    print(fermat_factorization(10000000000000000000000))
    end = time.time()
    print(f'Time: {end - start:.2f}s\n')
    x = "1.01"

