from euclicid import extendedeuclicid

def modular_linear_eq_solver(a, b, n):
    result = []
    [d, xp, yp] = extendedeuclicid(a, n)
    if b % d == 0:
        x0 = (xp * (b // d)) % n
        for i in range(0, d):
            result.append((x0 + i*(n // d)) % n)
    else:
        result = ["no solutions"]
    return result
print(modular_linear_eq_solver(14, 30, 100))
