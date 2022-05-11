from dataclasses import dataclass
from random import choice, seed

import time

class ElCurve:
    pass

class EcPoint:
    pass

class EcPointAffine:
    pass

class EcPointProjective:
    pass

@dataclass
class ElCurve:
    a:          int
    b:          int 
    order:      int
    field_size: int
    basepoint:  EcPoint

    def __init__(self, a, b, order, field_size, basepoint):
        self.a = a
        self.b = b
        self.order = order
        self.field_size = field_size
        basepoint.set_curve(self)
        self.basepoint = basepoint

        self.basepoint_affine = EcPointAffine(basepoint.x, basepoint.y, basepoint.z, basepoint.ec, False)
        self.basepoint_projective = EcPointProjective(basepoint.x, basepoint.y, basepoint.z, basepoint.ec)
@dataclass
class EcPoint:
    x:      int
    y:      int
    z:      int
    ec: ElCurve = None
  
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z and self.inf == other.inf

    def set_curve(self, ec):
        self.ec = ec

    def __repr__(self):
        return f'({self.x}, {self.y}, {self.z})'

@dataclass
class EcPointAffine(EcPoint):
    #if point is in infinity
    inf:    bool = False

    def __repr__(self):
        return f'({self.x}, {self.y}, {self.z}, {self.inf})'


    def __neg__(self):
        return EcPointAffine(self.x, -self.y, self.z, self.ec, self.inf) if not self.inf else self

    def __add__(self, other):
        return self._point_add(self, other)

    def __mul__(self, other: int):
        if self.ec:
            return fast_multiply(self, other, self.ec.a, self.ec.field_size)
        else:
            return fast_multiply(self, other, 0, None)

    def _point_add(self, P, Q):
        a = self.ec.a if self.ec else 0
        p = self.ec.field_size if self.ec else None

        #check if P is in infinity
        if (P.inf):
            return Q
        
        #check if q is in infinity
        if (Q.inf):
            return P

        #if they are equal, double
        if (Q == P):
            return self._point_double(P)

        #if p == -q return point in infinity
        if (P == -Q):
            return EcPointAffine(0, 1, 0, inf=True)

        alpha = (Q.y - P.y) * pow(Q.x - P.x, -1, p)
        beta = P.y - alpha * P.x 

        x = ((alpha ** 2 - (Q.x + P.x))) 
        y = ((- (beta + alpha * x)))
        if p != None:
            x = x % p
            y = y % p

        return EcPointAffine(x, y, self.z, self.ec, self.inf)


    def _point_double(self, P):
        a = self.ec.a if self.ec else 0
        p = self.ec.field_size if self.ec else None
        #return point at infinity
        if P.inf:
            return P

        alpha = ( 3 * P.x ** 2 + a) * pow(2 * P.y, -1, p) 
        beta = P.y - alpha * P.x

        x = ((alpha ** 2 - 2 * P.x) % p)
        y = ((- (beta + alpha * x)) % p)

        return EcPointAffine(x, y, self.z, self.ec, self.inf)

class EcPointProjective(EcPoint):
    @property
    def inf(self):
        return self.x == 0 and self.z == 0

    @inf.setter
    def inf(self):
        self.x = 0
        self.y = 1
        self.z = 0 

    def __neg__(self):
        return EcPointProjective(self.x, -self.y, self.z, self.ec) if not self.inf else self

    def __add__(self, other):
        return self._point_add(self, other)

    def __mul__(self, other: int):
        if self.ec:
            return fast_multiply(self, other, self.ec.a, self.ec.field_size)
        else:
            return fast_multiply(self, other, 0, None)

    def _point_add(self, P, Q):
        a = self.ec.a if self.ec else 0
        p = self.ec.field_size if self.ec else None

        #check if P is in infinity
        if (P.inf):
            return Q
        
        #check if q is in infinity
        if (Q.inf):
            return P

        #assign u1,u2,s1,s2
        u1 = P.x * Q.z
        u2 = Q.x * P.z
        s1 = P.y * Q.z
        s2 = Q.y * P.z

        #if they are equal, double
        if ( (s1 == s2) and (u1 == u2) ):
            return self._point_double(P)

        #if p == -q return point in infinity
        #if (P == -Q):
        #    return EcPointAffine(0, 1, 0, inf=True)

        #assign RR and PP
        PP = u2 - u1
        RR = s2 - s1

        #calculate alpha == RR/PP
        alpha = (RR * pow(PP, -1, p)) % p
        
        #assign w
        w = P.z * Q.z

        x = PP * (w * pow(RR, 2 ,p) - (u1 + u2) * pow(PP, 2, p))
        y = RR * (-2 * w * pow(RR, 2, p) +3 * (u1 + u2) * pow(PP,2,p)) - (s1 + s2) * pow(PP, 3, p)

        if y % 2 == 1:
            y += p
        
        y = y//2

        z = pow(PP, 3, p) * w

        x = x % p
        y = y % p
        z = z % p

        return EcPointProjective(x, y, z, self.ec)


    def _point_double(self, P):
        a = self.ec.a if self.ec else 0
        p = self.ec.field_size if self.ec else None
        #return point at infinity
        if P.inf:
            return P
        if P.y % p == 0 % p :
            return EcPointProjective(0, 1, 0, self.ec)

        #assign w and s
        w = ((3 * pow(P.x, 2, p) + a * pow(P.z, 2, p)) ) % p
        s = (2 * P.y * P.z) % p

        # alpha = w/s mod p
        alpha = (w * pow(s, -1, p)) % p

        # assign B and h
        B = (2 * s * P.x * P.y) % p
        h = (pow(w, 2, p) - 2 * B) % p

        x = (h * s) % p
        y = (w * (B-h) - 2 * pow(s*P.y, 2, p)) % p
        z = (pow(s, 3, p))

        return EcPointProjective(x, y, z, self.ec)


"""
double-and-add algorithm for fast multiplication
"""
def fast_multiply(P: EcPoint, s: int, a: int, p: int) -> EcPoint:
    b_bin = format(s, "b")
    
    R = P
    for bit in b_bin[1:]:
        R = R + R
        if bit == "1":
            R = R + P
    return R

EC = {
    6: ElCurve(34, 10, 41, 47, EcPoint(30,26,z=1)),
    20: choice([
        ElCurve(189977, 474810, 975619, 975259, EcPoint(889952, 89780, 1)),
    ]),
    30: choice([
        ElCurve(160392446, 570426436, 683656579, 683645161, EcPoint(268700899, 53018703, 1)),
    ]),
    40: choice([
    #    ElCurve(261420417203, 19624690186, 563542841207, 563542741247,  EcPoint(235493738933, 343108369384, 1)),
    #    ElCurve(518248484388, 284048620065, 743926924691, 743926365691,  EcPoint(208953459335, 177899353565, 1)),
    #    ElCurve(89469309901, 703706107702, 870032728657, 870031698089, EcPoint(550378620978, 19921237070, 1)),
        ElCurve(48739303714, 137819470075, 672538592221, 672539724617, EcPoint(305462345939, 467466940940, 1))
    ]),
    41: ElCurve(928679669726, 890471861893, 1116129833641, 1116129559639, EcPoint(821270767984, 195684596746, 1)),

    60: ElCurve(181423078584416465, 112783863405657942, 612356410521995059 ,612356409660197287, EcPoint(63852232286812509, 270435989459221281,1)),
}

if __name__ == "__main__":
    seed(time.time())
    ec = EC[6]

    ec.basepoint = EcPointAffine(ec.basepoint.x, ec.basepoint.y, ec.basepoint.z, ec.basepoint.ec, False)

    print(ec.basepoint)
    print("Eliptic curve: " + str(ec))
    print(ec.basepoint, -ec.basepoint)
    print(EcPointAffine(0, 123, 1, True), -EcPointAffine(0, 1230, 1, True))
    print(ec.basepoint + ec.basepoint + ec.basepoint)
    point = EcPointAffine(0, 1230, 1, ec, True)
    print(ec.basepoint + point)
    print(point + point)
    point = -ec.basepoint
    print(point)
    print(ec.basepoint + point)
    print(ec.basepoint * 2)
    print(ec.basepoint * 3)
    print((point + point) * 2)
    assert (point + point) + (point + point) == point * 4, "Assertion error, (P + P) + (P + P) != 4P"
    assert (point + point) + (point + point) == point * 2 + point * 2, "Assertion error, (P + P) + (P + P) != 2P + 2P"
    assert point * 2 + point * 3 == (point * 2) * 2 + point, "Assertion error, P*2 + P*3 != (P*2)*2 + 1"

    ec.basepoint = EcPointProjective(ec.basepoint.x, ec.basepoint.y, ec.basepoint.z, ec.basepoint.ec)

    print(ec.basepoint)
    print("Eliptic curve: " + str(ec))
    print(ec.basepoint, -ec.basepoint)
    print(ec.basepoint + ec.basepoint + ec.basepoint)
    point = EcPointProjective(0, 1230, 0, ec)
    point2 = EcPointProjective(1000, 0, 1000, ec)
    print(ec.basepoint + point)
    print(point + point)
    print(point.inf)
    print(point2.inf)
    print((point2 + point2))
    point = -ec.basepoint
    print(point)
    print(ec.basepoint + point)
    print(ec.basepoint * 2)
    print(ec.basepoint * 3)
    print((point + point) * 2)
    assert (point + point) + (point + point) == point * 4, "Assertion error, (P + P) + (P + P) != 4P"
    assert (point + point) + (point + point) == point * 2 + point * 2, "Assertion error, (P + P) + (P + P) != 2P + 2P"
    assert point * 2 + point * 3 == (point * 2) * 2 + point, "Assertion error, P*2 + P*3 != (P*2)*2 + 1"