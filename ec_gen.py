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
    6: ElCurve(34, 10, 41, 47, (30,26)),
    20: choice([
        ElCurve(189977, 474810, 975619, 975259, (889952, 89780)),
    ]),
    30: choice([
        ElCurve(160392446, 570426436, 683656579, 683645161, (268700899, 53018703)),
    ]),
    40: choice([
        ElCurve(261420417203, 19624690186, 563542841207, 563542741247,  (235493738933, 343108369384)),
        ElCurve(518248484388, 284048620065, 743926924691, 743926365691,  (208953459335, 177899353565)),
        ElCurve(89469309901, 703706107702, 870032728657, 870031698089, (550378620978, 19921237070)),
        ElCurve(48739303714, 137819470075, 672538592221, 672539724617, (305462345939, 467466940940))
    ]),
    41: ElCurve(928679669726, 890471861893, 1116129833641, 1116129559639, (821270767984, 195684596746)),

    60: ElCurve(181423078584416465, 112783863405657942, 612356410521995059 ,612356409660197287, (63852232286812509, 270435989459221281)),
}

def point_add(P, Q, a, p):
    if (P[0] == Q[0]):
        return point_double(P, a, p)
    
    alpha = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p)
    beta = P[1] - alpha * P[0] 

    x = ((alpha ** 2 - (Q[0] + P[0])) % p) 
    y = ((- (beta + alpha * x)) % p)

    return (x, y)


def point_double(P, a, p):
    if P[1] == 0:
        return P

    alpha = ( 3 * P[0] ** 2 + a) * pow(2 * P[1], -1, p) 
    beta = P[1] - alpha * P[0]

    x = ((alpha ** 2 - 2 * P[0]) % p)
    y = ((- (beta + alpha * x)) % p)

    return (x, y)

if __name__ == "__main__":
    seed(time.time())
    ec = EC[40]

    print("Eliptic curve: " + str(ec))