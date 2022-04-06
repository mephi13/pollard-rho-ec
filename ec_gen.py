from dataclasses import dataclass
from random import choice, seed
import time

@dataclass
class ElCurve:
    a:          int
    b:          int 
    order:      int
    field_size: int
    basepoint:  tuple

EC = {
    40: choice([
        ElCurve(261420417203, 19624690186, 563542841207, 563542741247,  (235493738933, 343108369384, 1)),
        ElCurve(518248484388, 284048620065, 743926924691, 743926365691,  (208953459335, 177899353565, 1)),
    ]),
    41: ElCurve(928679669726, 890471861893, 1116129833641, 1116129559639, (821270767984, 195684596746, 1)),
    #42: ElCurve(),
    #43: ElCurve(),

    60: ElCurve(181423078584416465, 112783863405657942, 612356410521995059 ,612356409660197287, (63852232286812509, 270435989459221281, 1)),
}

def point_addition(P, Q, a):
    if (P == Q)
        return point_doubling(P, a)
    
    alpha = (Q[1] - P[1])/(Q[0] - P[0])
    beta = P[1] - alpha * P[0] 

    x = pow(alpha, 2) - (Q[0] + P[0])
    y = - (beta + alpha * x)
    return (x, y)

def point_doubling(P, a):
    if p[1] == 0:
        return P

    alpha = ( 3 * P[0] ** 2 + a) / (2 * P[1]) 
    beta = P[1] - alpha * P[0]

    x = pow(alpha, 2) - 2 * P[0]
    y = - (beta + alpha * x)
    return (x, y)

if __name__ == "__main__":
    seed(time.time())
    ec = EC[40]

    print("Eliptic curve: " + str(ec))