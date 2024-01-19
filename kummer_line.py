"""
Implementation of the Kummer Line of Montgomery Curves and the corresponding 
Kummer Points used for x-only Montgomery curve arithmetic

===========================================================================

INFO: Construction


xP = K(X, Z)

Where x(P) = (X : Z) is the x-coordinate in projective XZ-coordinates

A KummerPoint can also be made straight from an elliptic curve point

E = EllipticCurve(F, [0,A,0,1,0])
P = E.random_point()

K = KummerLine(E)
xP = K(P)

===========================================================================

INFO: Usage

The methods of the KummerLine class are fairly straight-forward. Currently
missing is a check for whether two curves are isomorphic over the base field.

For this, we need an isomorphism between KummerLines which is a TODO.

For the points, scalar multiplication is performed by n*xP

Additionally, one can call `xP.double()` to perform x-only point addition
and xP.add(xQ, xPQ) to perform differential addition to recover xP + xQ
where xPQ = x(P - Q).

The 3 point ladder `xQ.ladder_3_pt(xP, xPQ, m) computes xP + [m]xQ

xP.multiples() generates values [l]xP by repeated differential addition. This
is used for isogeny computations where we want to collect the the first d points
for an isogeny of degree ell = 2d+1. 
"""

r"""
Kummer Line class for short Weierstrass elliptic curves

<Paragraph description>

EXAMPLES:

We construct a Kummer line over an already defined elliptic curve in short Weierstrass form::

    sage: E = EllipticCurve(GF(101), [2, 3])
    sage: KummerLine(E)
    Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101

We can also give the curve coefficients directly with a base ring::

    sage: KummerLine(QQ, [4, 5/6])
    Kummer line of the elliptic curve y^2 = x^3 + 4*x + 5/6 over Rational Field

TODO KummerPoint examples

AUTHORS:

- YOUR NAME (2005-01-03): initial version

- person (date in ISO year-month-day format): short desc

"""

# ****************************************************************************
#       Copyright (C) 2013 YOUR NAME <your email>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import cypari2

pari = cypari2.Pari()

from sage.all import cached_method, Integer, EllipticCurve

from sage.structure.element import RingElement
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field


class KummerLine:
    r"""
    Kummer line class for a short Weierstrass elliptic curve

    EXAMPLES::

    A KummerLine can be constructed straight from a short Weierstrass curve::
        sage: E = EllipticCurve(GF(101), [2, 3])
        sage: KummerLine(E)
        Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101

    Or, it can be constructed from the curve coefficients a and b, in which
    case a base ring must be specified::
        sage: KummerLine(QQ, [4, 5/6])
        Kummer line of the elliptic curve y^2 = x^3 + 4*x + 5/6 over Rational Field
    """

    def __init__(self, *args):
        r"""
        Construct a Kummer line from short Weierstrass coefficients or an elliptic curve
        shaped that way.

        INPUT::

            One of these two work:

            - ``E`` -- an EllipticCurve object in short Weierstrass form (y^2 = x^3 + a*x + b)

            OR

            - ``K`` -- a ring

            - ``[a, b]`` -- list or tuple of short Weierstrass curve coefficients

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: KummerLine(E)
            Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
            sage: KummerLine(QQ, [4, 5/6])
            Kummer line of the elliptic curve y^2 = x^3 + 4*x + 5/6 over Rational Field

        The base ring need not be a field::

            sage: KummerLine(IntegerModRing(15), [4, 7])
            Kummer line of the elliptic curve y^2 = x^3 + 4*x + 7 over Ring of integers modulo 15
        """
        self._curve = None

        # Allow the creation of the Kummer line from an EllipticCurve
        if len(args) == 1:
            (curve,) = args
            if not isinstance(curve, EllipticCurve_generic):
                raise TypeError("E is not an elliptic curve.")
            ainvs = curve.a_invariants()
            a, b = ainvs[3], ainvs[4]
            if ainvs != (0, 0, 0, a, b):
                raise ValueError("Must use short Weierstrass model.")
            self._curve = curve
            self._base_ring = curve.base_ring()

        # Allow the creation of the Kummer line from a base field and coeffs.
        elif len(args) == 2:
            base_ring, curve_constants = args
            # Extract curve constants
            if len(curve_constants) == 2:
                a, b = curve_constants
            else:
                raise ValueError(
                    "The coefficients must be a tuple [a, b] representing \
                    the two parameters of the elliptic curve in short Weierstrass form."
                )
            self._base_ring = base_ring
        else:
            raise ValueError(
                "A Kummer Line must be constructed from either a short Weierstrass curve, or\
                    a base field and tuple representing the coefficients [a, b]."
            )

        # init variables
        self._a = self._base_ring(a)
        self._b = self._base_ring(b)
        self._a, self._b = pari(self._a), pari(self._b)

        # Make sure the curve is not singular
        if self.discriminant() == 0:
            raise ValueError(
                f"Constants {curve_constants} do not define an elliptic curve in short \
                Weierstrass form."
            )

    def __eq__(self, other):
        r"""
        Test equality of Kummer lines by checking if the coefficients are the same.

        EXAMPLES::

            sage: KummerLine(QQ, [3, 4]) == KummerLine(QQ, [2, 3])
            False

            sage: KummerLine(QQ, [3, 4]) == KummerLine(QQ, [3, 4])
            True

        Base ring must be the same::
            sage: KummerLine(QQ, [3, 4]) == KummerLine(GF(7), [3, 4])
            False
        """

        if self.base_ring() != other.base_ring():
            return False
        return self._a == other._a and self._b == other._b

    def __repr__(self):
        r"""
        String representation of a Kummer line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.__repr__()
            'Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101'

            sage: K = KummerLine(QQ, [4, 5/6]); K.__repr__()
            'Kummer line of the elliptic curve y^2 = x^3 + 4*x + 5/6 over Rational Field'
        """

        return f"Kummer line of the elliptic curve y^2 = x^3 + {self._a.sage()}*x + {self._b.sage()} over {self.base_ring()}"

    def __call__(self, coords):
        r"""
        Create a point from the coordinates.

        INPUT::

            - ``coords`` - either a point P on EllipticCurve or a list or tuple (X, Z) where P = (X : * : Z); Z is optional

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)

        Point as a parameter::
            sage: K(P)
            Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101

        Same output if we just give the XZ-coordinates::
            sage: K([95, 1]) == K(P)
            True

        Z is optional::
            sage: K(95) == K(P)
            True

            sage: K([95]) == K(P)
            True
        """
        return KummerPoint(self, coords)

    def base_ring(self):
        r"""
        Return the base ring of the Kummer Line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: K.base_ring()
            Finite Field of size 101
        """
        return self._base_ring

    def extract_constants(self):
        r"""
        Return the curve coefficients a, b as a tuple.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: K.extract_constants()
            (2, 3)
        """
        # TODO should the example reflect that a and b are PARI elements?
        return self._a, self._b

    def zero(self):
        r"""
        Return the identity point on the Kummer Line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: K.zero()
            Kummer Point [1 : 0] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        """
        return self(None)

    def curve(self):
        r"""
        Lift the Kummer Line to an elliptic curve as a
        SageMath EllipticCurve.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: K.curve() == E
            True

            sage: K = KummerLine(QQ, [3, 4/5])
            sage: K.curve()
            Elliptic Curve defined by y^2 = x^3 + 3*x + 4/5 over Rational Field
        """

        if not self._curve:
            self._curve = self.short_weierstrass_curve()
        return self._curve

    @cached_method
    def short_weierstrass_curve(self):
        r"""
        Return the short Weierstrass curve associated with the Kummer line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: K.short_weierstrass_curve() == E
            True

            sage: K = KummerLine(QQ, [3, 4/5])
            sage: K.short_weierstreass_curve()
            Elliptic Curve defined by y^2 = x^3 + 3*x + 4/5 over Rational Field
        """

        F = self.base_ring()
        a, b = self.extract_constants()
        return EllipticCurve(F, [a, b])

    @cached_method
    def j_invariant(self):
        r"""
        Compute the j-invariant of the Kummer line using the formula `j(E) = 1728\frac{4a^3}{4a^3 + 27b^2}`.

        EXAMPLES::

            sage: E = EllipticCurve(j=42)
            sage: K = KummerLine(E)
            sage: K.j_invariant() == 42
            True

            sage: K = KummerLine(GF(101), [2, 3])
            sage: K.j_invariant()
            74
        """
        a, b = self.extract_constants()

        j_num = 4 * a**3
        j_den = j_num + 27 * b**2
        j_num = 1728 * j_num
        return j_num / j_den

    @cached_method
    def discriminant(self):
        r"""
        Compute the discriminant of the Kummer line using the formula `\Delta(E) = -16(4a^3 + 27b^2)`.

        EXAMPLES::

            sage: E = EllipticCurve(j=42)
            sage: K = KummerLine(E)
            sage: K.discriminant()
            -541067272574976

            sage: K = KummerLine(GF(101), [2, 3])
            sage: K.discriminant()
            44
        """
        a, b = self.extract_constants()
        return -16 * (4 * a**3 + 27 * b**2)

    def isogeny(self, S, P):
        XP, _ZP = P.XZ()
        a4, a6 = self.extract_constants()
        if self.zero() in S:
            S.remove(self.zero())
        v, w = 0, 0
        alpha = XP
        for Q in S:
            XQ, _ZQ = Q.XZ()
            gQx = 3 * XQ**2 + a4
            # if Q.double() == self.zero():
            if Q.double() == self.zero():
                # # print((2 * Q).x(), (Q.double().x()))
                # print(type(2 * Q), 2 * Q)
                # print(type(Q.double()), Q.double())
                # print(Q)
                vQ = gQx
            else:
                vQ = 2 * gQx
            uQ = 4 * (XQ**3 + a4 * XQ + a6)
            v += vQ
            w += uQ + XQ * vQ
            alpha += vQ / (XP - XQ) + uQ / (XP - XQ) ** 2
        A4 = a4 - 5 * v
        A6 = a6 - 7 * w
        K2 = KummerLine(self.base_ring(), [A4, A6])
        return K2, K2(alpha.sage())


class KummerPoint:
    r"""
    Kummer line point class

    EXAMPLES::

        sage: E = EllipticCurve(GF(101), [2, 3])
        sage: K = KummerLine(E)
        sage: P = E(95, 52)

    A KummerPoint can be constructed from coordinates::
        sage: xP = K([95, 1])
        Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101

    Or, it can be constructed from the curve coefficients a and b, in which
    case a base ring must be specified::
        sage: KummerLine(QQ, [4, 5/6])
        Kummer line of the elliptic curve y^2 = x^3 + 4*x + 5/6 over Rational Field
    """

    def __init__(self, parent, coords):
        r"""
        Create a point from the coordinates.

        INPUT::

            - ``coords`` - either a point P on EllipticCurve or a list or tuple (X, Z) where P = (X : * : Z); Z is optional

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)

        Point as a parameter::
            sage: K(P)
            Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101

        Same output if we just give the XZ-coordinates::
            sage: K([95, 1]) == K(P)
            True

        Z is optional::
            sage: K(95) == K(P)
            True

            sage: K([95]) == K(P)
            True
        """
        # Ensure the parent is the right type
        if not isinstance(parent, KummerLine):
            raise TypeError("not a Weierstrass Kummer line")

        R = parent.base_ring()

        # Point at infinity
        if coords is None:
            coords = (Integer(1), Integer(0))
        # Construct point from P on an elliptic curve in Montgomery form
        elif isinstance(coords, EllipticCurvePoint_field):
            # Make sure point's parent curve matches with Kummer Line
            a, b = parent._a, parent._b
            assert coords.curve().a_invariants() == (0, 0, 0, a, b)
            coords = coords[0], coords[2]
        # Construct from X coordinate only
        elif isinstance(coords, RingElement):
            coords = (coords,)
        # Construct from a tuple (X : Z)
        else:
            coords = tuple(coords)

        # Sanitise the input coordinates
        if len(coords) == 1:
            coords += (R.one(),)
        if len(coords) != 2:
            raise ValueError("not a point on ℙ¹")
        coords = tuple(map(pari, map(R, coords)))

        # TODO: we should make sure the coordinates
        #       are on the curve!
        # Write something like `parent.is_x_coord(X, Z)`
        # which checks if x^3 + A*x^2 + x is a square in
        # the base field...
        self._base_ring = R
        self._parent = parent
        self._X, self._Z = coords

    def __repr__(self):
        return f"Kummer Point [{self._X.sage()} : {self._Z.sage()}] on {self._parent}"

    def __bool__(self):
        """
        A point represents False if it is the point at infinity and
        True otherwise
        """
        return bool(self._Z)

    def __eq__(self, other):
        """
        Equality of two Kummer Points
        """
        if not isinstance(other, KummerPoint):
            raise ValueError("Can only compare equality between to Kummer Points")
        if self._parent != other._parent:
            return False
        return self._X * other._Z == other._X * self._Z

    def is_zero(self):
        """
        A Kummer Point is considered Zero if it is the identity point
        on the parent curve
        """
        return self._Z == 0

    def base_ring(self):
        """
        Get the base ring of the Kummer Point coordinates
        """
        return self._base_ring

    def parent(self):
        """
        Get the Kummer Line of which this point is constructed on
        """
        return self._parent

    def XZ(self):
        """
        Return the projective (X : Z) coordinates of the point
        """
        return self._X, self._Z

    def x(self):
        r""" """
        if not self._Z:
            raise ValueError("The identity point has no valid x-coordinate")
        if self._Z == 1:
            return self._base_ring(self._X)
        return self._base_ring(self._X / self._Z)

    @cached_method
    def curve_point(self):
        """
        Deterministically lift an x-coordinate
        taking the smallest y-coordinate as the
        chosen root.
        """
        # Get the Montgomery curve and constant A
        L = self.parent()
        E = L.curve()
        a, b = L._a, L._b

        # Compute y2, assume x is a valid coordinate
        x = self.x()
        y2 = x**3 + a * x + b
        y = y2.sqrt()
        return E(x, y)

    # =================================== #
    # Addition and multiplication helpers #
    # =================================== #

    @staticmethod
    def xDBL(X, Z, a, b2, b4):
        """
        function for Montgomery doubling with projective curve constant

        Input:  projective point P = (X:Z), curve constants (A:C)
        Output: projective point [2]P = (X2:Z2)

        Cost: 4M + 2S + 8a
        """

        XX = X * X
        ZZ = Z * Z
        t0 = X + Z
        t1 = t0 * t0 - XX - ZZ
        t1 = t1 + t1
        t0 = a * ZZ
        X2 = XX - t0
        X2 = X2 * X2 - b2 * t1 * ZZ
        Z2 = XX + t0
        ZZ = ZZ * ZZ
        Z2 = t1 * Z2 + b4 * ZZ

        return X2, Z2

    @staticmethod
    def xADD(XP, ZP, XQ, ZQ, xPQ, zPQ, a, b):
        """
        function for Montgomery differential addition

        Input:  projective coordinates P = (XP : ZP),
                Q=(XQ : ZQ), and their difference
                x(P-Q) = (xPQ : zPQ)
        Output: coordinates of sum P + Q = (XQP : ZQP)

        Cost: 4M + 2S + 6a
        """
        T1 = XP * XQ
        T2 = ZP * ZQ
        T3 = XP * ZQ
        T4 = ZP * XQ
        T5 = a * T2
        T6 = T1 - T5
        T7 = T6 * T6
        T8 = b * T2
        T9 = T8 + T8
        T9 = T9 + T9
        T10 = T3 + T4
        T11 = T9 * T10
        T12 = T7 - T11
        XQP = zPQ * T12
        T13 = T3 - T4
        T13 = T13 * T13
        ZQP = xPQ * T13

        return XQP, ZQP

    @staticmethod
    def xDBLADD(XP, ZP, XQ, ZQ, xPQ, zPQ, a, b4):
        """
        function for step in Montgomery ladder
        simultaneous doubling and differential addition

        Input: projective coordinates P=(XP:ZP) and Q=(XQ:ZQ),
               projective difference P-Q=(xPQ:zPQ) and
               curve constant A24/C24=(A+2C)/4C.
        Output: projective coordinates of 2P=(X2P:Z2P)
                and Q+P=(XQP:ZQP)

        Cost: 8M + 4S + 8A
        """
        XX = XP**2
        ZZ = ZP**2
        aZZ = a * ZZ
        t0 = XP + ZP
        t1 = t0**2
        t2 = t1 - XX
        E = t2 - ZZ
        t3 = XX - aZZ
        t4 = t3**2
        t5 = E * ZZ
        t6 = b4 * t5
        X2P = t4 - t6
        t7 = XX + aZZ
        t8 = ZZ**2
        t9 = b4 * t8
        t10 = E * t7
        t11 = t10 + t10
        Z2P = t11 + t9
        A = XP * XQ
        B = ZP * ZQ
        C = XP * ZQ
        D = XQ * ZP
        t12 = a * B
        t13 = A - t12
        t14 = C + D
        t15 = t13**2
        t16 = B * t14
        t17 = b4 * t16
        t18 = t15 - t17
        XQP = zPQ * t18
        t19 = C - D
        t20 = t19**2
        ZQP = xPQ * t20

        return X2P, Z2P, XQP, ZQP

    # =================================== #
    # Addition and multiplication methods #
    # =================================== #
    def _double(self):
        """
        Returns [2] self
        """
        X, Z = self.XZ()
        a, b = self._parent.extract_constants()
        b2 = b + b
        b4 = b2 + b2
        X2, Z2 = self.xDBL(X, Z, a, b2, b4)
        return self._parent((X2, Z2))

    def _double_iter(self, n):
        """
        Returns [2^k] self
        """
        X, Z = self.XZ()
        a, b = self._parent.extract_constants()
        b2 = b + b
        b4 = b2 + b2
        for _ in range(n):
            X, Z = self.xDBL(X, Z, a, b2, b4)
        return self._parent((X, Z))

    def double(self):
        """
        Wrapper function which deals with the doubling of
        the identity

        Returns [2] * self
        """
        # Deal with identity
        if not self._Z:
            return self
        return self._double()

    def double_iter(self, n):
        """
        Wrapper function which deals with the repeated
        doubling

        Returns [2^n] * self

        This avoids the ADD part of xDBLADD, and so is
        faster when we know our scalar is a power of two
        """
        # Deal with identity
        if not self._Z:
            return self
        return self._double_iter(n)

    def _add(self, Q, PQ):
        """
        Performs differential addition assuming
        P, Q and PQ are all not the point at
        infinity
        """
        XP, ZP = self.XZ()
        XQ, ZQ = Q.XZ()
        XPQ, ZPQ = PQ.XZ()

        a, b = self._parent.extract_constants()

        X_new, Z_new = self.xADD(XP, ZP, XQ, ZQ, XPQ, ZPQ, a, b)
        return self._parent((X_new, Z_new))

    def add(self, Q, PQ):
        """
        Function to perform differential addition and
        handle the cases when P, Q or PQ are the points
        at infinity
        """
        # Adding O + Q = Q
        if not self._Z:
            return Q

        # Adding P + O = P
        if not Q._Z:
            return self

        # Difference is the identity
        # so P = Q and P+Q = [2]P
        if not PQ._Z:
            return self._double()

        return self._add(Q, PQ)

    def __mul__(self, m):
        """
        Montgomery-ladder to compute [m]P

        Input: coordinates of P=(XP:ZP)
               scalar factor m, curve constants (A:C)
        Output: KummerPoint [m]P=(X0:Z0)
        """
        if not isinstance(m, (int, Integer)):
            try:
                m = Integer(m)
            except:
                raise TypeError(f"Cannot coerce input scalar {m = } to an integer")

        # If m is zero, return identity
        if not m:
            return self.parent().zero()

        # [m]P = [-m]P for x-only
        m = abs(m)

        # Extract base field and coefficients
        R = self.base_ring()
        XP, ZP = self.XZ()

        # Initialise for loop
        X0, Z0 = R.one(), R.zero()
        X1, Z1 = XP, ZP

        # Converting parameters for projective DBLADD -> (A24:C24)=(A+2C:4C)
        a, b = self.parent().extract_constants()
        b4 = b + b
        b4 = b4 + b4
        a, b4 = pari(a), pari(b4)
        # A24 = C + C
        # C24 = pari(A24 + A24)
        # A24 = pari(A24 + A)

        # Montgomery-ladder
        for bit in bin(m)[2:]:
            if bit == "0":
                X0, Z0, X1, Z1 = self.xDBLADD(X0, Z0, X1, Z1, XP, ZP, a, b4)
            else:
                X1, Z1, X0, Z0 = self.xDBLADD(X1, Z1, X0, Z0, XP, ZP, a, b4)

        return self._parent((X0, Z0))

    def __rmul__(self, m):
        return self * m

    def __imul__(self, m):
        self = self * m
        return self

    def ladder_3_pt(self, xP, xPQ, m):
        """
        Function to compute x(P + [m]Q) using x-only
        arithmetic. Very similar to the Montgomery ladder above

        Note: self = xQ
        """
        if not isinstance(m, (int, Integer)):
            try:
                m = Integer(m)
            except:
                raise TypeError(f"Cannot coerce input scalar {m = } to an integer")

        # If m is zero, return xP
        if not m:
            return xP

        # [m]P = [-m]P for x-only
        m = abs(m)

        # Converting parameters for projective DBLADD -> (A24:C24)=(A+2C:4C)
        a, b = self.parent().extract_constants()
        b4 = b + b
        b4 = b4 + b4
        a, b4 = pari(a), pari(b4)

        # Extract out coordinates
        XQ, ZQ = self.XZ()
        XP, ZP = xP.XZ()
        XPQ, ZPQ = xPQ.XZ()

        # Montgomery-ladder
        for bit in bin(m)[:1:-1]:
            if bit == "1":
                XQ, ZQ, XP, ZP = self.xDBLADD(XQ, ZQ, XP, ZP, XPQ, ZPQ, a, b4)
            else:
                XQ, ZQ, XPQ, ZPQ = self.xDBLADD(XQ, ZQ, XPQ, ZPQ, XP, ZP, a, b4)
        return self._parent((XP, ZP))

    def multiples(self):
        """
        A generator of points [l]P for self = P
        Stops when it has generated the full subgroup generated by P
        (without the identity point).

        NOTE: this is implemented to make Vélu-like computations easy
        """
        yield self
        R = self.double()
        # Order 2 case
        if not R:
            return

        # Odd order case
        Q = self
        while R:
            yield R
            S = R.add(self, Q)
            Q, R = R, S

        return
