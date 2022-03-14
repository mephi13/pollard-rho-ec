from sympy import randprime, isprime, pollard_rho as sympy_pollard_rho, nextprime
from random import randrange
import time


"""   
Generate a strong prime p between min and max, such that p = 2p_t + 2 where p, p_t are prime
"""
def gen_strong_prime(min, max, randomize_bin_len = False):
    p = 0
    min_v = min
    max_v = max
    while not isprime(p):
        if randomize_bin_len:
            min_v = randrange(min,max)
            max_v = min_v + 1
        p_t = randprime(2 ** min_v, 2 ** max_v)
        p = 2*p_t + 1
    return (p, p_t)

"""   
Generate base of the DLP g_t and g. Order of g_t must be equal to p_t
"""
def gen_g(p, p_t = None):
    g = 0
    g_t = 1
    while g_t == 1:
        g = randrange(0, p)
        g_t = pow(g,2,p) 
    if p_t:
        #Order of g_t should be p_t    
        assert(pow(g_t, p_t, p) == 1)
    return (g, g_t)

"""
generate y in <g> such that y = g^x mod p
x will be then calculated using pollard rho alghorhitm 
"""
def gen_y(g_t, p_t, p):
    x = randrange(0, p_t)
    y = pow(g_t, x, p)
    return (x, y)

"""
Calculate next parameters Xi+1, ai+1 and bi+1 based on the mapping function fi(Xi, ai, bi) -> (Xi+1, ai+1, bi+1)
"""
def f_mapping(Xi, ai, bi, g_t, p, y, p_t):
    case = Xi % 3

    #Xi belongs to S0
    if case == 1: 
        Xj = g_t * Xi % p
        aj = ai + 1 % p_t
        bj = bi 
    #Xi belongs to S1
    elif case == 2:
        Xj = pow(Xi, 2, p)
        aj = 2*ai % p_t
        bj = 2*bi % p_t
    #Xi belongs to S2
    else:
        Xj = Xi *y % p
        aj = ai
        bj = bi + 1 % p_t
    return (Xj, aj, bj)

def pollard_rho(g_t, p, p_t, y):
    i = 0
    #Slow "Tortoise"
    T, a, b = (1, 0, 0)
    
    while True:
        #Fast "Hare"
        H, g, d = (T, a, b)

        i += 1
        T, a, b = f_mapping(T, a, b, g_t, p, y, p_t)
        H, g, d = f_mapping(T, a, b, g_t, p, y, p_t) 
        while (T != H % p):
            i += 1
            T, a, b = f_mapping(T, a, b, g_t, p, y, p_t)
            H, g, d = f_mapping(H, g, d, g_t, p, y, p_t)
            H, g, d = f_mapping(H, g, d, g_t, p, y, p_t)
            
        #Now T = H, so a + xb = g + xd mod p_t
        if b != d % p_t:
            x = (a-g)*pow(d - b, -1, p_t) % p_t
            return (x, i)
        else:
            a = randrange(0, p_t)
            b = randrange(0, p_t)
            T = pow(g_t, a + b, p)
            print(f"Starting pollard-rho again, with alpha = {a}, beta = {b}, T = {T}")

def main():
    min = 40
    max = 59
    
    p, p_t = gen_strong_prime(min, max, True)
    print(f'p is a {p.bit_length()} bit number')
    print("p = %s, p_t = %s" % (p, p_t))

    g, g_t = gen_g(p, p_t)
    print("g = %s, g_t = %s" % (g , g_t))
    print(f"Check order: g_t^p_t mod p = {(pow(g_t, p_t, p))}" )

    real_x, y = gen_y(g_t, p_t, p)
    print(f"Real x = {real_x}, \ny = g^x mod p = {y}")
    print(f'Estimated maximum iterations order of magnitude sqrt(p_t) = {round(pow(p_t, 0.5))}')

    start_time = time.time()
    calc_x, i = pollard_rho(g_t, p, p_t, y)
    calc_time = time.time() - start_time
    print(f'Calcutaion time: {calc_time}')
    print(f"Iterations: {i}")
    print(f"Real x = {real_x}, calculated x = {calc_x}")
    try:
        assert(pow(g_t, calc_x, p) == y)
        print("Real x and calculated x are the same")
    except:
        print("Assertion failure, calculated x does not equal real x")

if __name__=="__main__":
    main()