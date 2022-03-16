#!/home/mephi/.virtualenvs/pollard/bin/python
#for prime generation
from sympy import randprime, isprime, nextprime
#misc
import argparse
from random import randrange
import time


"""   
Generate a strong prime p between min and max, such that p = 2p_t + 2 where p, p_t are prime
"""
def gen_strong_prime(min, max, randomize_bin_len = False, p =None, p_t = None):
    if p:
        p_t = (p - 1) // 2
    elif p_t:
        p =  2*p_t + 1
    else:
        p = 0
        min_v = min - 2
        max_v = max - 2
        while not isprime(p):
            if randomize_bin_len:
                min_v = randrange(min - 2,max - 2)
                max_v = min_v + 1
            p_t = randprime(2 ** min_v, 2 ** max_v)
            p = 2*p_t + 1
    #print(f"Bit length: {p.bit_length()}")
    return (p, p_t)

"""   
Generate base of the DL g_t and g. Order of g_t must be equal to p_t
"""
def gen_g(p, g = None):
    if g != None:
        g_t = pow(g, 2, p) 
    else:
        g = 0
        g_t = 1
        while g_t == 1:
            g = randrange(1, p)
            g_t = pow(g,2,p) 
    return (g, g_t)

"""
generate y in <g> such that y = g^x mod p
x will be then calculated using pollard rho alghorhitm 
"""
def gen_y(g_t, p_t, p, real_x = None):
    x = randrange(0, p_t) if not real_x else real_x
    y = pow(g_t, x, p)
    return (x, y)

"""
Pollard-rho alghorhitm for DL problems
"""
def pollard_rho(g_t, p, p_t, y):
    def f_mapping(Xi, ai, bi):
        case = Xi % 3
        #Xi belongs to S0
        if case == 1: 
            Xj = (g_t * Xi) % p
            aj = (ai + 1) % p_t
            bj = bi 
        #Xi belongs to S1
        elif case == 2:
            Xj = pow(Xi, 2, p)
            aj = (2*ai) % p_t
            bj = (2*bi) % p_t
        #Xi belongs to S2
        else:
            Xj = (Xi *y) % p
            aj = ai
            bj = (bi + 1) % p_t
        return (Xj, aj, bj)

    i = 0
    #Slow "Tortoise"
    T, a, b = (1, 0, 0)
    #Fast "Hare"
    H, g, d = (1, 0, 0)

    while True:
        i += 1

        T, a, b = f_mapping(T, a, b)
        H, g, d = f_mapping(*f_mapping(H, g, d))

        if (T == H % p): 
            #Now T = H, so a + xb = g + xd mod p_t
            if  b != d % p_t:
                x = (a-g)*pow(d - b, -1, p_t) % p_t
                return (x, i)
            else:
                a = randrange(0, p_t)
                b = randrange(0, p_t)
                T = (pow(g_t, a, p) * pow(y, b, p)) % p
                H, g, d = (T, a, b)
                print(f"Algorithm unsuccessful - d == b\nStarting pollard-rho again with alpha = {a}, beta = {b}, T = {T}")

"""
Genearate an instance of an DL problem
"""
def generate_DLP_instance(min_p, max_p, p, p_t, g, g_t , real_x, y):
    try:
        if not (p and p_t):
            p, p_t = gen_strong_prime(min_p, max_p, True, p, p_t)
        assert (isprime(p) and isprime(p_t) and p == 2*p_t + 1), "p or p_t are not primes or p != 2*p_t +1, use different values"
        assert (p_t != 2), "p_t can't be equal to 2, use a different value"

        if not (g_t):
            g, g_t = gen_g(p, g=g)
        assert (pow(g_t, p_t, p) == 1), "Order of g_t should be equal to p_t, use different value"

        if not (y):
            real_x, y = gen_y(g_t, p_t, p, real_x)
        elif real_x:
            assert (pow(g_t, real_x, p) == y), "g_t^x mod p != y, use correct x value"

        return (p, p_t, g, g_t, real_x, y)

    except AssertionError as ex:
        raise Exception(f"Generating DLP failed, reason: \n{ex.args[0]}")

def main(min = None, max = None, p= None, p_t = None, g = None, g_t = None, real_x = None, y = None):
    try:
        nmin = 40 if not min else min
        max = 61 if not (max or min) else min + 1 if max==None else max + 1
        min = nmin

        p, p_t, g, g_t, real_x, y = generate_DLP_instance(min, max, p, p_t, g, g_t, real_x, y)

        print(f'p is a {p.bit_length()} bit number')
        print("p = %s, p_t = %s" % (p, p_t))
        print("g = %s, g_t = %s" % (g , g_t)) if g else print("g_t = %s" % (g_t)) 
        print(f"Check order of g_t: {g_t}^{p_t} mod {p} = {(pow(g_t, p_t, p))}" )
        print(f"Real x = {real_x}, \ny = {g_t}^{real_x} mod {p} = {y}") if real_x else print(f"y = {g_t}^x mod {p} = {y}")
        print(f'Estimated maximum iterations order of magnitude sqrt(p_t) = {round(pow(p_t, 0.5)):,}')

        start_time = time.time()
        calc_x, i = pollard_rho(g_t, p, p_t, y)
        calc_time = time.time() - start_time
        
        print(f'Calcutaion time: {round(calc_time, 1)} seconds')
        print(f"Iterations: {i:,}")
        print(f"Real x = {real_x}, calculated x = {calc_x}") if real_x else print(f"Calculated x = {calc_x}")

        if pow(g_t, calc_x, p) == y:
            print(f"Calculations correct, {g_t}^{calc_x} mod {p} = {y}")
        else:
            print(f"Calculations failed,  {g_t}^{calc_x} mod {p} != {y}. y might not be in <g>")
    except Exception as ex:
        print(ex.args[0])

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Polland-rho alghorhitm for DL problems in finite fields.')
    parser.add_argument("--min", type=int)
    parser.add_argument("--max", type=int)
    parser.add_argument("-p", type=int)
    parser.add_argument("-pt", type=int)
    parser.add_argument("-g", type=int)
    parser.add_argument("-gt", type=int)
    parser.add_argument("-x", type=int)
    parser.add_argument("-y", type=int)
    a = parser.parse_args()

    main(a.min, a.max, a.p, a.pt, a.g, a.gt, a.x, a.y)