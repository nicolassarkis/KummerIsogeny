from kummer_line import KummerLine
from sage.all import EllipticCurve, GF, randint
from benchmark_utils import compare_isogeny

# p = 65537
# F = GF(p)
#
# E = EllipticCurve(F, [3, 4])
# P = E.random_point()
# Q = E.random_point()
# K = KummerLine(E)
# KP = K(P)
# KQ = K(Q)
# KPQ = K(P - Q)
#
# n = 3481
# k = 13
# R = KP.double_iter(k)
# print((2**k * P)[0] == R.x())
#
# # K2 = KummerLine_montgomery(F, [3, 4])
# # print(K2.a())
# # print(K2.a()._coeff_repr())
#
# print(K(None))
# print(E.j_invariant() == K.j_invariant())
#
# R = n * KP
# print((n * P)[0] == R.x())
#
# print("-------------")
#
# F = GF(101)
# E = EllipticCurve(F, [2, 3])
# P = E.random_point()
# K = KummerLine(E)
# xP = K([P[0], P[2]])
# print(P)
# print(xP)

p = 65537
F = GF(p)

A = F.random_element()
B = F.random_element()
E = EllipticCurve(F, [A, B])
K = KummerLine(E)
N = E.order()
k = N.divisors()[1:][randint(0, len(N.divisors()) - 2)]

P0 = E.gens()[0]
P = N / k * P0
assert P.order() == k
# while True:
#     # l = random_prime(37)
#     l = 8 * 11 * 107
#     if N % l == 0:
#         P = E.random_element()
#         P = N / l * P
#         assert l * P == E(0) and P != E(0)
#         print("l = ", l)
#         S = []
#         i = 1
#         while i * P != E(0):
#             S = S + [i * P]
#             i += 1
#         S = S[: ceil(len(S) // 2)]
#         print(S)
#         KS = [K(P) for P in S]
#         break

Q = E.random_point()
xP = K(P)
xQ = K(Q)

compare_isogeny(P, Q, xP, xQ, k)

# K2, KR = K.iso_velu(KS, KQ)
# XKR, ZKR = KR.XZ()
# phi = E.isogeny(S)
# E2 = phi.codomain()
# R = phi(Q)
# xR, yR = R.xy()

# print(E2.j_invariant() == K2.j_invariant())
# print(xR * ZKR == XKR)
