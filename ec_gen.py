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

@dataclass
class EcPoint:
    x:      int
    y:      int
    #if point is in infinity
    inf:    bool = False
    def __repr__(self):
        return f'({self.x}, {self.y})'

    def __neg__(self):
        return EcPoint(self.x, -self.y, False) if not self.inf else self


EC = {
    6: ElCurve(34, 10, 41, 47, EcPoint(30,26)),
    20: choice([
        ElCurve(189977, 474810, 975619, 975259, EcPoint(889952, 89780)),
    ]),
    30: choice([
        ElCurve(160392446, 570426436, 683656579, 683645161, EcPoint(268700899, 53018703)),
    ]),
    40: choice([
        ElCurve(261420417203, 19624690186, 563542841207, 563542741247,  EcPoint(235493738933, 343108369384)),
        ElCurve(518248484388, 284048620065, 743926924691, 743926365691,  EcPoint(208953459335, 177899353565)),
        ElCurve(89469309901, 703706107702, 870032728657, 870031698089, EcPoint(550378620978, 19921237070)),
        ElCurve(48739303714, 137819470075, 672538592221, 672539724617, EcPoint(305462345939, 467466940940))
    ]),
    41: ElCurve(928679669726, 890471861893, 1116129833641, 1116129559639, EcPoint(821270767984, 195684596746)),

    60: ElCurve(181423078584416465, 112783863405657942, 612356410521995059 ,612356409660197287, EcPoint(63852232286812509, 270435989459221281)),
}

def point_add(P, Q, a, p):
    #check if P is in infinity
    if (P.inf):
        return Q
    
    #check if q is in infinity
    if (Q.inf):
        return P

    #if they are equal, double
    if (Q == P):
        return point_double(P, a, p)

    #if p == -q return point in infinity
    if (P == -Q):
        return EcPoint(0, 0, inf=True)

    alpha = (Q.y - P.y) * pow(Q.x - P.x, -1, p)
    beta = P.y - alpha * P.x 

    x = ((alpha ** 2 - (Q.x + P.x)) % p) 
    y = ((- (beta + alpha * x)) % p)

    return EcPoint(x, y)


def point_double(P, a, p):
    #return point at infinity
    if P.inf:
        return P

    alpha = ( 3 * P.x ** 2 + a) * pow(2 * P.y, -1, p) 
    beta = P.y - alpha * P.x

    x = ((alpha ** 2 - 2 * P.x) % p)
    y = ((- (beta + alpha * x)) % p)

    return EcPoint(x, y)


if __name__ == "__main__":
    seed(time.time())
    ec = EC[40]

    print("Eliptic curve: " + str(ec))
    print(ec.basepoint, -ec.basepoint)
    print(EcPoint(0, 123, True), -EcPoint(0, 1230, True))