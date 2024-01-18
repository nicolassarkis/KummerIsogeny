from kummer_line import KummerLine_weierstrass
from sage.all import EllipticCurve, GF

p = 65537
F = GF(p)

E = EllipticCurve(F, [3, 4])
P = E.random_point()
Q = E.random_point()
K = KummerLine_weierstrass(E)
KP = K(P)
KQ = K(Q)
KPQ = K(P - Q)

n = 3481
k = 13
R = KP.double_iter(k)
print((2**k * P)[0] == R.x())

# K2 = KummerLine_montgomery(F, [3, 4])
# print(K2.a())
# print(K2.a()._coeff_repr())

print(K(None))
print(E.j_invariant() == K.j_invariant())

R = n * KP
print((n * P)[0] == R.x())
