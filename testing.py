from kummer_line import KummerLine  # noqa
from sage.all import EllipticCurve, GF, randint, Integer  # noqa
from benchmark_utils import compare_isogeny  # noqa


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

# c = 0
# flag = True
#
# while c < 2000 and flag:
#     A = F.random_element()
#     B = F.random_element()
#     E = EllipticCurve(F, [A, B])
#     N = E.order()
#
#     P0 = E.gens()[0]
#     while P0.order() == N:
#         A = F.random_element()
#         B = F.random_element()
#         E = EllipticCurve(F, [A, B])
#         N = E.order()
#         P0 = E.gens()[0]
#
#     K = KummerLine(E)
#     k = P0.order().divisors()[1:][randint(0, len(P0.order().divisors()) - 2)]
#     P = P0.order() / k * P0
#     assert P.order() == k
#
#     ker = [E(0)]
#     for _ in range(P.order() - 1):
#         ker += [ker[-1] + P]
#
#     Q = E.random_point()
#     while Q in ker:
#         Q = E.random_point()
#     xP = K(P)
#     xQ = K(Q)
#
#     flag = compare_isogeny(P, Q, xP, xQ, k)
#     c += 1
#     if c % 100 == 0:
#         print(c)

# A = F.random_element()
# B = F.random_element()
# E = EllipticCurve(F, [A, B])
# if E.order() % 2 == 0:
#     print("allo")
#     P = E.gens()[0]
#     P = P.order() / 2 * P
#     assert P.order() == 2
# else:
#     P = E.random_point()
# Q = E.random_point()
# K = KummerLine(E)
# xP, xQ, xPQ = K(P), K(Q), K(P - Q)
# print(P, xP.x())
# print(Q, xQ.x())
# print(P - Q, xPQ.x())
#
# dbl = xP.xDBL(xP.x(), 1, A, 2 * B, 4 * B)
# print(dbl[1] == 0)
#
# add = xP.xADD(xP.x(), 1, xQ.x(), 1, xPQ.x(), 1, A, B)
# print((P + Q) == add[0] / add[1])
#
# dbladd = xP.xDBLADD(xP.x(), 1, xQ.x(), 1, xPQ.x(), 1, A, 4 * B)
# print(dbladd[1] == 0, dbladd[2] / dbladd[3] == add[0] / add[1])
#
# print((2 * xP).is_zero())
#
# print(P)
# print(2 * P)
# print(xP)
# print(2 * xP)
# print(xP.double())

E = EllipticCurve(F, [60860, 10913])
K = KummerLine(E)
P = E.lift_x(Integer(0))
xP = K(P)
print(xP)
print("aa", xP.xDBLADD(0, 1, 1, 0, 0, 1, 60860, 43652))
print(2 * xP)
print(xP.double().x() == (2 * P)[0])
