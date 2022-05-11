#!/bin/python
#for prime generation
from sympy import randprime, isprime, nextprime
from ec_gen import EC, EcPointAffine, ElCurve, EcPoint
#misc
import argparse
from termcolor import colored
from random import randrange
import time


def gen_Y_ec(P, q, real_x = None):
    s = randrange(2, q - 1) if not real_x else real_x
    Y = P * s
    return (s, Y)

"""
Pollard-rho algorithm for DL problems
"""
def pollard_rho(Y, P, q):
    def f_mapping(Xi, ai, bi):
        case = Xi.x % 3
        #Xi belongs to S0
        if case == 1: 
            Xj = Xi + P
            aj = (ai + 1) % q
            bj = bi 
        #Xi belongs to S1
        elif case == 2:
            Xj = Xi + Xi
            aj = (2*ai) % q
            bj = (2*bi) % q
        #Xi belongs to S2
        else:
            Xj = Xi + Y
            aj = ai
            bj = (bi + 1) % q
        return (Xj, aj, bj)

    i = 0
    #Slow "Tortoise"
    T, alpha, beta = (P, 1, 0)
    #Fast "Hare"
    H, gamma, delta = (P, 1, 0)

    while True:
        i += 1

        T, alpha, beta = f_mapping(T, alpha, beta)
        H, gamma, delta = f_mapping(*f_mapping(H, gamma, delta))

        #print(f"T: {T}, H: {H}")
        if (T == H): 
            #Now T = H, so aP + bQ = gP + dQ = -(a-g)/(d-b)P = Q
            if  beta != delta % q:
                #print(alpha, gamma, delta, beta)
                s = (alpha-gamma)*pow(delta - beta, -1, q) % q
                return (s, i)
            else:
                alpha = randrange(0, q)
                beta = randrange(0, q)
                T = (P * alpha) + (Y * beta) 
                H, gamma, delta = (T, alpha, beta)
                print(f"Algorithm unsuccessful - d == b\nStarting pollard-rho again with alpha = {alpha}, beta = {beta}, T = {T}")

"""
Genearate an instance of an DL problem
"""
def generate_DLP_instance(n_bits, real_s, Y):
    try:
        #Get pregenerated EC
        assert n_bits in EC, f"There is no pre-generated {n_bits} bits long EC, use a different value"
        ec = EC[n_bits]
        #Generate Y
        if not (Y):
            real_s, Y = gen_Y_ec(ec.basepoint_affine, ec.order, real_s)
        elif real_s:
            assert (P * real_s == Y), "P*s mod p != y, use correct x value"

        return (ec, real_s, Y)

    except Exception as ex:
        raise Exception(f"Generating DLP failed, reason: \n{ex.args[0]}")


def main(n_bits = None, real_s = None, Y = None):
    try:
        #init params
        n_bits = 40 if not n_bits else n_bits
        ec, real_s, Y = generate_DLP_instance(n_bits, real_s, Y)

        print(f'p is a {ec.field_size.bit_length()} bit number')
        print(f'EC: {str(ec)}')
        print(f"Real s = {colored(real_s, 'magenta')},") 
        print(f"Y = {colored(ec.basepoint, 'green')} * {colored(real_s, 'magenta')}"\
            f" mod {colored(ec.field_size, 'yellow')}\nY = {colored(Y, 'blue')}")\
            if real_s else print(f"Y = {colored(ec.basepoint, 'green')}*s mod {colored(ec.field_size, 'yellow')}\nY = {colored(Y, 'blue')}")
        print(f'Estimated maximum iterations order of magnitude sqrt(q) = {round(pow(ec.order, 0.5)):,}')


        #Pollard-rho calculations + time checking
        print(f"Starting pollard-rho(Y, P, q, p, a, b) = ")
        print(f"Pollard-Rho({colored(Y, 'blue')}, {colored(ec.basepoint, 'green')}, {colored(ec.order, 'cyan')}, {colored(ec.field_size, 'yellow')}, {colored(ec.a, 'white')}, {colored(ec.b, 'white')})")
        start_time = time.time()
        calc_s, i = pollard_rho(Y, ec.basepoint_affine, ec.order)
        calc_time = time.time() - start_time
        print(f'Calculation time: {round(calc_time, 1)} seconds')
        print(f"Iterations: {i:,}")
        print(f"Real s = {colored(real_s, 'magenta')}, calculated s = {colored(calc_s, 'red')}") if real_s else print(f"Calculated s = {colored(calc_s, 'red')}")
    

        #Check if real s and calc s are the same
        if (ec.basepoint_affine * calc_s) == Y:
            print(f"{colored(ec.basepoint, 'green')}*{colored(calc_s, 'red')} mod {colored(ec.field_size, 'yellow')} = {colored(Y, 'blue')}")
            print(f"{colored('Calculations correct!', 'green')}")
        else:
            print(f"{colored(ec.basepoint, 'green')}*{colored(calc_s, 'red')} mod {colored(ec.field_size, 'yellow')} = {colored(Y, 'blue')}")
            print(f"{colored('Calculations failed!', 'red')}")
    except Exception as ex:
        print(ex.args[0])

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Polland-rho alghorhitm for DL problems over Eliptic Curve.')
    parser.add_argument("-n", type=int)
    parser.add_argument("-y", nargs='+',type=int)
    parser.add_argument("-s", type=int)
    a = parser.parse_args()
    if a.y != None:
        assert len(a.y) >= 2, "Supply at least 2 coordinates for Y"
        y = EcPointAffine(a.y[0], a.y[1], 1) if len(a.y) <3 else EcPoint(a.y[0], a.y[1], a.y[2]) 
    else:
        y = a.y

    #print(fast_multiply(7, 8))
    #print(fast_multiply(EC[40].basepoint, 16, EC[40].a ))

    main(a.n, a.s, y)