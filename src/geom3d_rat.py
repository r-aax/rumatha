"""
3D geometry in rational coordinates.
"""

import numpy as np
from fractions import Fraction as Fr
import geom2d_rat
import matplotlib.pyplot as plt

#===================================================================================================

class Point:
    """
    Point with three coordinates.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, x=Fr(0), y=Fr(0), z=Fr(0)):
        """
        Constructor.

        Parameters
        ----------
        x : Fraction
            x-coordinate.
        y : Fraction
            y-coordinate.
        z : Fraction
            z-coordinate.
        """

        self.x = x
        self.y = y
        self.z = z

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_real_coords(x, y, z, denom):
        """
        Constructor from real coordinates.

        Parameters
        ----------
        x : float
            Real coordinate X.
        y : float
            Real coordinate Y.
        z : float
            Real coordinate Z.
        denom : int
            Denominator.

        Returns
        -------
        Point
            Constructed point.
        """

        return Point(Fr(int(x * denom), denom),
                     Fr(int(y * denom), denom),
                     Fr(int(z * denom), denom))

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_real_array(ar, denom):
        """
        Constructor from real array.

        Parameters
        ----------
        ar : np.array
            Array of coordinates.
        denom : int
            Denominator.

        Returns
        -------
        Point
            Constructed point.
        """

        return Point.from_real_coords(ar[0], ar[1], ar[2], denom)

    #-----------------------------------------------------------------------------------------------

    def get_real_array(self):
        """
        Get real array from rat coordinates.

        Returns
        -------
        np.array
            Array of coordinates.
        """

        coords = [self.x, self.y, self.z]

        return [f.numerator / f.denominator for f in coords]

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, p):
        """
        Check equal with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if equal to another point,
            False - otherwise.
        """

        return (self.x == p.x) and (self.y == p.y) and (self.z == p.z)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, p):
        """
        Check not equal to another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if not equal to another point.
            False - otherwise.
        """

        return not (self == p)

    #-----------------------------------------------------------------------------------------------

    def __ge__(self, p):
        """
        Check for GE with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if self >= p,
            False - otherwise.
        """

        if self.x > p.x:
            return True
        elif self.x < p.x:
            return False
        elif self.y > p.y:
            return True
        elif self.y < p.y:
            return False
        else:
            return self.z >= p.z

    #-----------------------------------------------------------------------------------------------

    def __gt__(self, p):
        """
        Check for GT with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if self > p,
            False - otherwise.
        """

        if self.x > p.x:
            return True
        elif self.x < p.x:
            return False
        elif self.y > p.y:
            return True
        elif self.y < p.y:
            return False
        else:
            return self.z > p.z

    #-----------------------------------------------------------------------------------------------

    def __le__(self, p):
        """
        Check for LE with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if self <= p,
            False - otherwise.
        """

        return not (self > p)

    #-----------------------------------------------------------------------------------------------

    def __lt__(self, p):
        """
        Check for LT with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - it self < p,
            False - otherwise.
        """

        return not (self >= p)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'P({self.x}, {self.y}, {self.z})'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', size=20):
        """
        Draw on plit.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        size : int
            Size.
        """

        # Ignore no size points.
        if size > 0:

            #Ignore z coordinate.
            plt.scatter(self.x, self.y, color=color, s=size)

    #-----------------------------------------------------------------------------------------------

    def __sub__(self, p):
        """
        Two points difference.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        Vector
            Two points difference.
        """

        return Vector(self.x - p.x, self.y - p.y, self.z - p.z)

    #-----------------------------------------------------------------------------------------------

    def is_between(self, A, B):
        """
        Check if point is between two points.

        Parameters
        ----------
        A : Point
            Point.
        B : Point
            Point.

        Returns
        -------
        bool
            True - if point is between A and B,
            False - otherwise.
        """

        return ((A <= self) and (self <= B)) or ((B <= self) and (self <= A))

    #-----------------------------------------------------------------------------------------------

    def projection_OXY(self):
        """
        Projection on OXY (ignore Z coord).

        Returns
        -------
        geom2d_rat.Point
            Projection.
        """

        return geom2d_rat.Point(self.x, self.y)

    #-----------------------------------------------------------------------------------------------

    def projection_OXZ(self):
        """
        Projection on OXZ (ignore Y coord).

        Returns
        -------
        geom2d_rat.Point
            Projection.
        """

        return geom2d_rat.Point(self.x, self.z)

    #-----------------------------------------------------------------------------------------------

    def projection_OYZ(self):
        """
        Projection of OYZ (ignore X coord).

        Returns
        -------
        geom2d_rat.Point
            Projection.
        """

        return geom2d_rat.Point(self.y, self.z)

#===================================================================================================

class Vector(Point):
    """
    Class vector.
    """

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'V({self.x}, {self.y}, {self.z})'

    #-----------------------------------------------------------------------------------------------

    def __add__(self, v):
        """
        Add two vectors.

        Parameters
        ----------
        v : Vector.
            Vector.

        Returns
        -------
        Vector
            New vector.
        """

        return Vector(self.x + v.x, self.y + v.y, self.z + v.z)

    #-----------------------------------------------------------------------------------------------

    def __sub__(self, v):
        """
        Subtract two vectors.

        Parameters
        ----------
        v : Vector
            Vector.

        Returns
        -------
        Vector
            New vector.
        """

        return Vector(self.x - v.x, self.y - v.y, self.z - v.z)

    #-----------------------------------------------------------------------------------------------

    def is_null(self):
        """
        Check if vector is null.

        Returns
        -------
        bool
            True - if vector is null,
            False - otherwise.
        """

        return (self.x == 0) and (self.y == 0) and (self.z == 0)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def dot(v1, v2):
        """
        Dot product.

        Parameters
        ----------
        v1 : Vector
            First vector.
        v2 : Vector
            Second vector.

        Returns
        -------
        Fraction
            Result.
        """

        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def vector_product(a, b):
        """
        Vector product.

        Parameters
        ----------
        a : Vector
            First vector.
        b : Vector
            Second vector.

        Returns
        -------
        Vector
            Result vector.
        """

        return Vector(a.y * b.z - a.z * b.y,
                      a.z * b.x - a.x * b.z,
                      a.x * b.y - a.y * b.x)

#===================================================================================================

class Points:
    """
    Points.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.items = []

    #-----------------------------------------------------------------------------------------------

    def __getitem__(self, i):
        """
        Get i-th point.

        Parameters
        ----------
        i : int
            Index.

        Returns
        -------
        Point
            Point.
        """

        return self.items[i]

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Points{self.items}'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        size : int
            Size.
        """

        for p in self.items:
            p.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def count(self):
        """
        Count of points.

        Returns
        -------
        int
            Count of points.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    def add(self, p):
        """
        Add new point.

        Parameters
        ----------
        p : Point
            Point.
        """

        self.items.append(p)

    #-----------------------------------------------------------------------------------------------

    def add_unique(self, p):
        """
        Add new unique point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point was added,
            False - if point was not added.
        """

        for pi in self.items:
            if p == pi:
                return False

        self.add(p)

        return True

    #-----------------------------------------------------------------------------------------------

    def sort(self):
        """
        Sort points.
        """

        self.items.sort()

#===================================================================================================

class Line:
    """
    Line in space.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, x0, y0, z0, m, n, p):
        """
        Line constructor.

        x = x0 + tm
        y = y0 + tn
        z = z0 + tp

        Parameters
        ----------
        x0 : Fraction
            X coordinate of base point.
        y0 : Fraction
            Y coordinate of base point.
        z0 : Fraction
            Z coordinate of base point.
        m : Fraction
            Parameter for X direction.
        n : Fraction
            Parameter for Y direction.
        p : Fraction
            Parameter for Z direction.
        """

        # Normalize vector.
        if m != 0:
            n = n / m
            p = p / m
            m = Fr(1)
        elif n != 0:
            p = p / n
            n = Fr(1)
        else:
            assert p != 0
            p = Fr(1)

        # Normalize point.
        if m != 0:
            t = -x0 / m
            y0 = y0 + t * n
            z0 = z0 + t * p
            x0 = Fr(0)
        elif n != 0:
            t = -y0 / n
            x0 = x0 + t * m
            z0 = z0 + t * p
            y0 = Fr(0)
        else:
            assert p != 0
            t = -z0 / p
            x0 = x0 + t * m
            y0 = y0 + t * n
            z0 = Fr(0)

        self.P0 = Point(x0, y0, z0)
        self.v = Vector(m, n, p)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_point_and_vector(p, v):
        """
        Constructor from point and vector.

        Parameters
        ----------
        p : Point
            Base point.
        v : Vector
            Direction vector.

        Returns
        -------
        Line
            Constructed line.
        """

        return Line(p.x, p.y, p.z, v.x, v.y, v.z)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_points(p1, p2):
        """
        Constructor from two points.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point
            Second point.

        Returns
        -------
        Line
            Constructed line.
        """

        return Line.from_point_and_vector(p1, p2 - p1)

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, line):
        """
        Check equal.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if lines are equal,
            False - otherwise.
        """

        return (self.P0 == line.P0) and (self.v == line.v)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, line):
        """
        Check not equal.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if lines are not equal,
            False - otherwise.
        """

        return not (self == line)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Line(x = {self.x0} + t * {self.m}, '\
               f'y = {self.y0} + t * {self.n}, '\
               f'z = {self.z0} + t * {self.p})'

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if line has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if line has point,
            False - otherwise.
        """

        x, y, z = p.x, p.y, p.z
        x0, y0, z0, m, n, p = self.P0.x, self.P0.y, self.P0.z, self.v.x, self.v.y, self.v.z

        if m != 0:
            t = (x - x0) / m
            return (y == y0 + t * n) and (z == z0 + t * p)
        elif n != 0:
            t = (y - y0) / n
            return (x == x0 + t * m) and (z == z0 + t * p)
        elif p != 0:
            t = (z - z0) / p
            return (x == x0 + t * m) and (y == y0 + t * n)
        else:
            assert False

#===================================================================================================

class Segment:
    """
    Segment in space.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, A, B):
        """
        Constructor by points.

        Parameters
        ----------
        A : Point
            A Point.
        B : Point
            B Point.
        """

        # Construction of zero length segment is forbidden.
        assert A != B

        self.A = A
        self.B = B

        # Keep point in sorted way.
        self.sort_points()

        # Construct line for this segment.
        self.line = Line.from_points(self.A, self.B)

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, s):
        """
        Check equal to another segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if equal to another segment,
            False - otherwise.
        """

        # Actually we need only (self.A == s.A) and (self.B == s.B)
        # since points are sorted.
        return ((self.A == s.A) and (self.B == s.B)) or ((self.A == s.B) and (self.B == s.A))

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, s):
        """
        Check not equal to another segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if not equal to another segment,
            False - otherwise.
        """

        return not (self == s)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Segm[{self.A}, {self.B}]'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Line width.
        size : int
            Size.
        """

        # We ignore z coordinate, draw only on OXY plane.
        x = [self.A.x, self.B.x]
        y = [self.A.y, self.B.y]

        # Do not draw segment with zero width.
        if linewidth != '0':
            plt.plot(x, y, color=color, linewidth=linewidth)

        # Draw ends.
        self.A.draw(plt, color=color, size=size)
        self.B.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def mod2(self):
        """
        Square of module.

        Returns
        -------
        Fraction
            Square of module.
        """

        return (self.A.x - self.B.x)**2 + (self.A.y - self.B.y)**2 + (self.A.z - self.B.z)**2

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if segment has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if segment has point.
            False - otherwise.
        """

        return self.line.is_have_point(p) and p.is_between(self.A, self.B)

    #-----------------------------------------------------------------------------------------------

    def sort_points(self):
        """
        Sort points.
        """

        if self.A > self.B:
            self.A, self.B = self.B, self.A

#===================================================================================================

class Segments:
    """
    Segments.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.items = []

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Line width.
        size : int
            Size.
        """

        for s in self.items:
            s.draw(plt, color=color, linewidth=linewidth, size=size)

    #-----------------------------------------------------------------------------------------------

    def count(self):
        """
        Count of segments.

        Returns
        -------
        int
            Count of segments.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    def add(self, s):
        """
        Add segment.

        Parameters
        ----------
        s : Segment
            Segment.
        """

        self.items.append(s)

    #-----------------------------------------------------------------------------------------------

    def add_unique(self, s):
        """
        Add unique segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if new segment was added,
            False - if segment was not added.
        """

        for si in self.items:
            if si == s:
                return False

        self.add(s)

        return True

    #-----------------------------------------------------------------------------------------------

    def points(self):
        """
        Get all points.

        Returns
        -------
        Points
            Points.
        """

        ps = Points()

        for s in self.items:
            ps.add_unique(s.A)
            ps.add_unique(s.B)

        return ps

    #-----------------------------------------------------------------------------------------------

    def sort(self, fun):
        """
        Sort set of segments.

        Parameters
        ----------
        fun : fun
            Key function.
        """

        self.items.sort(key=fun)

#===================================================================================================

class Plane:
    """
    Plane.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, a, b, c, d):
        """
        Constructor.

        Parameters
        ----------
        a : Fraction
            a coefficient.
        b : Fraction
            b coefficient.
        c : Fraction
            c coefficient.
        d : Fraction
            d coefficient.
        """

        self.a = a
        self.b = b
        self.c = c
        self.d = d

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_points(A, B, C):
        """
        Create plane from three points.

        Plane types.
        1) if a != 0 then a = 1
           x + by + cz + d = 0
        2) if a = 0 and b != 0 then b = 1
           b + cz + d = 0
        3) if a = 0 and b = 0 and c != 0 then c = 1
           z + d = 0
        4) if a = 0 and b = 0 and c = 0 then there is no plane.

        Parameters
        ----------
        A : Point
            A Point.
        B : Point
            B Point.
        C : Point
            C Point.

        Returns
        -------
        Plane
            Result plane.
        """

        x1, y1, z1, x2, y2, z2, x3, y3, z3 = A.x, A.y, A.z, B.x, B.y, B.z, C.x, C.y, C.z

        # src: https://guimc.bmstu.ru/wp-content/uploads/2018/11/lecture_2.1.pdf
        # vectors M1M  = (x - x1, y - y1, z - z1)
        #         M1M2 = (x2 - x1, y2 - y1, z2 - z1)
        #         M1M3 = (x3 - x1, y3 - y1, z3 - z1)
        # These vectors lie in one plane when they are complanar:
        #   | x - x1    y - y1    z - z1  |
        #   | x2 - x1   y2 - y1   z2 - z1 | = 0
        #   | x3 - x1   y3 - y1   z3 - z1 |
        x21, y21, z21 = x2 - x1, y2 - y1, z2 - z1
        x31, y31, z31 = x3 - x1, y3 - y1, z3 - z1

        #   | x - x1    y - y1    z - z1 |
        #   |  x21       y21       z21   | = 0
        #   |  x31       y31       z31   |
        # Calculate determinant:
        # (x - x1) * |y21 z21| - (y - y1) * |x21 z21| + (z - z1) * |x21 y21| = 0
        #            |y31 z31|              |x31 z31|              |x31 y31|
        dx = y21 * z31 - y31 * z21
        dy = -(x21 * z31 - x31 * z21)
        dz = (x21 * y31 - x31 * y21)

        # (x - x1) * dx + (y - y1) * dy + (z - z1) * dz = 0
        # x dx - x1 dx + y dy - y1 dy + z dz - z1 dz = 0
        # x dx + y dy + z dz + (-(x1 dx + y1 dy + z1 dz)) = 0
        a = dx
        b = dy
        c = dz
        d = -(x1 * dx + y1 * dy + z1 * dz)

        if a != 0:
            b = b / a
            c = c / a
            d = d / a
            a = Fr(1)
        elif b != 0:
            c = c / b
            d = d / b
            b = Fr(1)
        elif c != 0:
            d = d / c
            c = Fr(1)
        else:
            raise Exception(f'Plane can not be constructed from points {A}, {B}, {C}.')

        return Plane(a, b, c, d)

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, p):
        """
        Check equal with another plane.

        Parameters
        ----------
        p : Plane
            Plane.

        Returns
        -------
        bool
            True - if equal to another plane,
            False - otherwise.
        """

        return (self.a == p.a) and (self.b == p.b) and (self.c == p.c) and (self.d == p.d)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, p):
        """
        Check not equal with another plane.

        Parameters
        ----------
        p : Plane
            Plane.

        Returns
        -------
        bool
            True - if not equal to another plane,
            False - otherwise.
        """

        return not (self == p)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Plane({self.a} x + {self.b} y + {self.c} z + {self.d})'

    #-----------------------------------------------------------------------------------------------

    def val(self, p):
        """
        Value of point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        Fraction
            Value.
        """

        return self.a * p.x + self.b * p.y + self.c * p.z + self.d

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if plane has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if plane has point,
            False - otherwise.
        """

        return self.val(p) == 0

    #-----------------------------------------------------------------------------------------------

    def is_two_points_strong_on_one_side(self, p1, p2):
        """
        Check if two points strong on one side.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point
            Second point.

        Returns
        -------
        bool
            True - if two points on one side of plane,
            False - otherwise.
        """

        return self.val(p1) * self.val(p2) > 0

    #-----------------------------------------------------------------------------------------------

    def normal(self):
        """
        Get normal vector.

        Returns
        -------
        Vector
            Normal vector.
        """

        return Vector(self.a, self.b, self.c)

    #-----------------------------------------------------------------------------------------------

    def is_perpendicular_with_plane(self, pl):
        """
        Check if plane is perpendicular with plane.

        Parameters
        ----------
        pl : Plane
            Plane.

        Returns
        -------
        bool
            True - if perpendicular with plane,
            False - otherwise.
        """

        return Vector.dot(self.normal(), pl.normal()) == 0

    #-----------------------------------------------------------------------------------------------

    def is_intersects_with_segment(self, s):
        """
        Check if plane intersects with segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if plane intersects segment,
            False - otherwise.
        """

        return not self.is_two_points_strong_on_one_side(s.A, s.B)

    #-----------------------------------------------------------------------------------------------

    def intersection_with_segment(self, s):
        """
        Find plane and segment intersection.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        None
            If there is no intersection.
        Point
            If there is only point of intersection.
        Segment
            If segment lay in plane.
        """

        # Check for no intersection.
        if not self.is_intersects_with_segment(s):
            return None

        # Check if both ends of segment lie in plane.
        if self.is_have_point(s.A) and self.is_have_point(s.B):
            return s

        # Intersection point is inside of segment.
        return self.intersection_with_line(Line.from_points(s.A, s.B))

    #-----------------------------------------------------------------------------------------------

    def is_intersects_with_line(self, line):
        """
        Check if plane intersects with line.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if plane intersects with line,
            False - otherwise.
        """

        # see intersection_with_line

        return not (self.intersection_with_line(line) is None)

    #-----------------------------------------------------------------------------------------------

    def intersection_with_line(self, line):
        """
        Find intersection with line.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        None
            If there is no intersection.
        Point
            If plane and line intersects by point.
        Line
            If line lays in plane.
        """

        a, b, c, d = self.a, self.b, self.c, self.d
        x0, y0, z0, m, n, p = line.P0.x, line.P0.y, line.P0.z, line.v.x, line.v.y, line.v.z

        # ax + by + cz + d = 0
        #   x = x0 + tm
        #   y = y0 + tn
        #   z = z0 + tp

        # One of m, n, p is not zero.

        if m != 0:

            # case 1. m != 0
            #   t = (x - x0) / m, y = y0 + (n / m)(x - x0), z = z0 + (p / m)(x - x0)
            #   ax + b (y0 + (n / m)(x - x0)) + c(z0 + (p / m)(x - x0)) + d = 0
            #   ax + b y0 + (bn / m)(x - x0) + c z0 + (cp / m)(x - x0) + d = 0
            #   ax + b y0 + (bn / m)x - (bn / m)x0 + c z0 + (cp / m)x - (cp / m)x0 + d = 0
            #   x (a + (bn + cp) / m) + b y0 + c z0 + d - ((bn + cp) / m) x0 = 0
            #   q = (bn + cp) / m
            #   x (a + q) + b y0 + c z0 + d - q x0 = 0
            #   x = (q x0 - b y0 - c z0 - d) / (a + q)
            q = (b * n + c * p) / m
            up, dn = q * x0 - b * y0 - c * z0 - d, a + q

            if dn == 0:
                if up == 0:
                    return line
                else:
                    return None

            x = up / dn
            t = (x - x0) / m
            y = y0 + t * n
            z = z0 + t * p

            return Point(x, y, z)

        elif n != 0:

            # case 2. n != 0
            #   t = (y - y0) / n, x = x0 + (m / n)(y - y0), z = z0 + (p / n)(y - y0)
            #   a (x0 + (m / n)(y - y0)) + by + c (z0 + (p / n)(y - y0)) + d = 0
            #   a x0 + (am / n)(y - y0) + by + c z0 + (cp / n)(y - y0) + d = 0
            #   a x0 + (am / n)y - (am / n)y0 + by + c z0 + (cp / n)y - (cp / n)y0 + d = 0
            #   y ( b + (am + cp) / n) + a x0 + c z0 + d - ((am + cp) / n) y0 = 0
            #   q = (am + cp) / n
            #   y (b + q) + a x0 + c z0 + d - q y0 = 0
            #   y = (q y0 - a x0 - c z0 - d) / (b + q)
            q = (a * m + c * p) / n
            up, dn = q * y0 - a * x0 - c * z0 - d, b + q

            if dn == 0:
                if up == 0:
                    return line
                else:
                    return None

            y = up / dn
            t = (y - y0) / n
            x = x0 + t * m
            z = z0 + t * p

            return Point(x, y, z)

        else:

            assert p != 0

            # case 3. p != 0
            #   t = (z - z0) / p, x = x0 + (m / p)(z - z0), y = y0 + (n / p)(z - z0)
            #   a (x0 + (m / p)(z - z0)) + b (y0 + (n / p)(z - z0)) + cz + d = 0
            #   a x0 + (am / p)(z - z0) + b y0 + (bn / p)(z - z0) + cz + d = 0
            #   a x0 + (am / p)z + (am / p)z0 + b y0 + (bn / p)z + (bn / p)z0 + cz + d = 0
            #   z (c + (am + bn) / p) + a x0 + b y 0 + d - ((am + bn) / p) z0 = 0
            #   q = (am + bn) / p
            #   z (c + q) + a x0 + b y0 + d - q z0 = 0
            #   z = (q z0 - a x0 - b y0 - d) / (c + q)
            q = (a * m + b * n) / p
            up, dn = q * z0 - a * x0 - b * y0 - d, c + q

            if dn == 0:
                if up == 0:
                    return line
                else:
                    return None

            z = up / dn
            t = (z - z0) / p
            x = x0 + t * m
            y = y0 + t * n

            return Point(x, y, z)

#===================================================================================================

class Triangle:
    """
    Triangle in space.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, A, B, C):
        """
        Constructor.

        Parameters
        ----------
        A : Point
            A Point.
        B : Point
            B Point.
        C : Point
            C Point.
        """

        self.A = A
        self.B = B
        self.C = C
        self.points = [self.A, self.B, self.C]

        # Init sides at a time.
        self.AB = Segment(self.A, self.B)
        self.BC = Segment(self.B, self.C)
        self.AC = Segment(self.A, self.C)
        self.sides = [self.AB, self.BC, self.AC]

        # Plane.
        self.plane = Plane.from_points(self.A, self.B, self.C)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Tri<{self.A}, {self.B}, {self.C}>'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Color.
        size : int
            Size.
        """

        # First draw sides without points.
        for s in self.sides:
            s.draw(plt, color=color, linewidth=linewidth, size=0)

        # Then draw points over sides.
        for p in self.points:
            p.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def is_perpendicular_with_plane(self, pl):
        """
        Check if perpendicular with plane.

        Parameters
        ----------
        pl : Plane
            Plane.

        Returns
        -------
        bool
            True - if perpendicular with plane,
            False - otherwise.
        """

        return self.plane.is_perpendicular_with_plane(pl)

    #-----------------------------------------------------------------------------------------------

    def projection_OXY(self):
        """
        Projection on OXY (ignore Z coordinate).

        Returns
        -------
        geom2d_rat.Triangle
            Projection.
        """

        return geom2d_rat.Triangle(self.A.projection_OXY(),
                                   self.B.projection_OXY(),
                                   self.C.projection_OXY())

    #-----------------------------------------------------------------------------------------------

    def projection_OXZ(self):
        """
        Projection on OXZ (ignore Y coordinate).

        Returns
        -------
        geom2d_rat.Triangle
            Projection.
        """

        return geom2d_rat.Triangle(self.A.projection_OXZ(),
                                   self.B.projection_OXZ(),
                                   self.C.projection_OXZ())

    #-----------------------------------------------------------------------------------------------

    def projection_OYZ(self):
        """
        Projection on OYZ (ignore X coordinate).

        Returns
        -------
        geom2d_rat.Triangle
            Projection.
        """

        return geom2d_rat.Triangle(self.A.projection_OYZ(),
                                   self.B.projection_OYZ(),
                                   self.C.projection_OYZ())

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if triangle has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if triangle has point,
            False - otherwise.
        """

        if not self.is_perpendicular_with_plane(OXY):
            return self.projection_OXY().is_have_point(p.projection_OXY())
        elif not self.is_perpendicular_with_plane(OXZ):
            return self.projection_OXZ().is_have_point(p.projection_OXZ())
        elif not self.is_perpendicular_with_plane(OYZ):
            return self.projection_OYZ().is_have_point(p.projection_OYZ())
        else:
            assert False

    #-----------------------------------------------------------------------------------------------

    def intersection_with_segment(self, s):
        """
        Find intersection of triangle with segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        None
            There is non intersection.
        Point
            If there is only one point intersection.
        Segment
            If intersection is performed by segment.
        """

        # Find intersection of plane and segment.
        r = self.plane.intersection_with_segment(s)

        # No intersection segment with plane.
        if r is None:
            return None

        # There is intersection of plane and segment.
        if isinstance(r, Point):
            p = r
            if self.is_have_point(p):
                return p
            else:
                return None

        assert isinstance(r, Segment)

        # Find intersections of all triangle sides with segment.
        r1 = Intersection.segment_segment(self.AB, s)
        r2 = Intersection.segment_segment(self.BC, s)
        r3 = Intersection.segment_segment(self.AC, s)

        # If one of intersections r1, r2, r3 is segment,
        # then this segment is result of triangle and segment intersection.
        if isinstance(r1, Segment):
            return r1
        if isinstance(r2, Segment):
            return r2
        if isinstance(r3, Segment):
            return r3

        # Then we have to collect all points.
        ps = Points()
        if isinstance(r1, Point):
            ps.add_unique(r1)
        if isinstance(r2, Point):
            ps.add_unique(r2)
        if isinstance(r3, Point):
            ps.add_unique(r3)

        # Now result is None, Point or Segment.
        cnt = ps.count()
        if cnt == 0:
            return None
        elif cnt == 1:
            return ps[0]
        else:
            assert cnt == 2
            return Segment(ps[0], ps[1])

    #-----------------------------------------------------------------------------------------------

    def intersection_with_triangle(self, t):
        """
        Find intersection with another triangle.

        Parameters
        ----------
        t : Triangle
            Triangle.

        Returns
        -------
        None
            There is no intersection.
        Point
            There is one point intersection.
        Segment
            There is two points intersection.
        [Point]
            Set of points, which form convex full.
            May contain points from 0 to 6.
            3 - triangle,
            4 - convex quadrangle,
            5 - convex pentagon,
            6 - convex hexagon.
        """

        ps = Points()

        for s in self.sides:
            r = t.intersection_with_segment(s)
            if isinstance(r, Point):
                ps.add_unique(r)
            elif isinstance(r, Segment):
                ps.add_unique(r.A)
                ps.add_unique(r.B)

        for s in t.sides:
            r = self.intersection_with_segment(s)
            if isinstance(r, Point):
                ps.add_unique(r)
            elif isinstance(r, Segment):
                ps.add_unique(r.A)
                ps.add_unique(r.B)

        cnt = ps.count()

        if cnt == 0:
            return None
        elif cnt == 1:
            return ps[0]
        elif cnt == 2:
            return Segment(ps[0], ps[1])
        else:
            return ps

#===================================================================================================

class Intersection:
    """
    Intersection of two geometrical objects.
    """

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_point(p1, p2):
        """
        Intersection of two points.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point.
            Secod point.

        Returns
        -------
        None
            No intersection.
        Point
            Equal points.
        """

        # Points may ne equal or not.
        if p1 == p2:
            return p1
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_line(p, ln):
        """
        Intersection of point and line.

        Parameters
        ----------
        p : Point
            Point.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on line.
        """

        # Point is on line or there is no intersection.
        if ln.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_segment(p, s):
        """
        Intersection of point and segment.

        Parameters
        ----------
        p : Point
            Point.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on segment.
        """

        # Point is on segment or there is no intersection.
        if s.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_plane(p, pl):
        """
        Intersection of point and plane.

        Parameters
        ----------
        p : Point
            Point.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on plane.
        """

        # Point is on plane or there is no intersection.
        if pl.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_triangle(p, t):
        """
        Intersection of point and triangle.

        Parameters
        ----------
        p : Point
            Point.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Point is inside triangle.
        """

        # Point is in triangle or there is no intersection.
        if t.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_point(ln, p):
        """
        Intersection of line and point.

        Parameters
        ----------
        ln : Line
            Line.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on line.
        """

        return Intersection.point_line(p, ln)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_line(ln1, ln2):
        """
        Intersection of two lines.

        Parameters
        ----------
        ln1 : Line
            First line.
        ln2 : Line
            Second line.

        Returns
        -------
        None
            No intersection.
        Point
            Single points intersection.
        Line
            The same line.
        """

        # First check for equal.
        if ln1 == ln2:
            return ln1

        # Check for parallel lines.
        if Vector.vector_product(ln1.v, ln2.v).is_null():
            return None

        # Extract coefficients.
        x1, y1, z1 = ln1.P0.x, ln1.P0.y, ln1.P0.z
        m1, n1, p1 = ln1.v.x, ln1.v.y, ln1.v.z
        x2, y2, z2 = ln2.P0.x, ln2.P0.y, ln2.P0.z
        m2, n2, p2 = ln2.v.x, ln2.v.y, ln2.v.z

        # Linear equations system.
        # x1 + t1 * m1 = x2 + t2 * m2
        # y1 + t1 * n1 = y2 + t2 * n2
        # z1 + t1 * p1 = z2 + t2 * p2
        # Move all members to left.
        # t1 * m1 - t2 * m2 + (x1 - x2) = 0
        # t1 * n1 - t2 * n2 + (y1 - y2) = 0
        # t1 * p1 - t2 * p2 + (z1 - z2) = 0
        m2, n2, p2 = -m2, -n2, -p2
        dx, dy, dz = x1 - x2, y1 - y2, z1 - z2

        # System in simple form.
        # m1 * t1 + m2 * t2 + dx = 0 // (1)
        # n1 * t1 + n2 * t2 + dy = 0 // (2)
        # p1 * t1 + p2 * t2 + dz = 0 // (3)

        # Function for solving each system of 3.
        def slv2(m1, m2, dx, n1, n2, dy):
            # (1) * n1 - (2) * m1
            #   (m2 * n1 - n2 * m1) * t2 + (dx * n1 - dy * m1) = 0
            #   t2 = (dy * m1 - dx * n1) / q
            # (1) * n2 - (2) * m2
            #   (m1 * n2 - n1 * m2) * t1 + (dx * n2 - dy * m2) = 0
            #   t1 = (dx * n2 - dy * m2) / q
            q = m2 * n1 - n2 * m1
            if q == 0:
                return None
            else:
                return  (dx * n2 - dy * m2) / q, (dy * m1 - dx * n1) / q

        # Try to solve system of equations.
        r = slv2(m1, m2, dx, n1, n2, dy)
        if r is None:
            r = slv2(m1, m2, dx, p1, p2, dz)
            if r is None:
                r = slv2(n1, n2, dy, p1, p2, dz)
        t1, t2 = r

        # Try all equations.
        if (m1 * t1 + m2 * t2 + dx == 0) \
            and (n1 * t1 + n2 * t2 + dy == 0) \
            and (p1 * t1 + p2 * t2 + dz == 0):
            return Point(x1 + t1 * m1, y1 + t1 * n1, z1 + t1 * p1)
        else:
            return None

    #--------------------------------------------------------------------------------------

    @staticmethod
    def line_segment(ln, s):
        """
        Intersection of line and segment.

        Parameters
        ----------
        ln : Line
            Line.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment is on line.
        """

        # Find intersection of two lines.
        r = Intersection.line_line(ln, s.line)

        # No intersection.
        if r is None:
            return None

        # Whole segment.
        if isinstance(r, Line):
            return s

        # Intersection of two lines is point.
        assert isinstance(r, Point)

        p = r

        # If point on segment then this point is intersection.
        if s.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_plane(ln, pl):
        """
        Intersection of line and plane.

        Parameters
        ----------
        ln : Line
            Line.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Line
            Line is in plane.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_triangle(ln, t):
        """
        Intersection of line and triangle.

        Parameters
        ----------
        ln : Line
            Line.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_point(s, p):
        """
        Intersection of segment and point.

        Parameters
        ----------
        s : Segment
            Segment.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on segment.
        """

        return Intersection.point_segment(p, s)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_line(s, ln):
        """
        Intersection of segment and line.

        Parameters
        ----------
        s : Segment
            Segment.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment is on line.
        """

        return Intersection.line_segment(ln, s)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_segment(s1, s2):
        """
        Intersection of two segments.

        Parameters
        ----------
        s1 : Segment
            Segment.
        s2 : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        # Find intersection of two containing lines.
        r = Intersection.line_line(s1.line, s2.line)

        # No intersection.
        if r is None:
            return None

        # One point intersection.
        # Both segments must have it.
        # We know that point is placed on both lines.
        # So we have to check point for between ends check.
        if isinstance(r, Point):
            p = r
            if p.is_between(s1.A, s1.B) and p.is_between(s2.A, s2.B):
                return p
            else:
                return None

        # So both segments lie in one line.
        assert isinstance(r, Line)

        # Get ends.
        A1, B1 = s1.A, s1.B
        A2, B2 = s2.A, s2.B
        assert A1 <= B1
        assert A2 <= B2

        # No intersection check.
        # A1       B1       A2       B2
        # *--------*........*--------*
        # A2       B2       A1       B1
        #*---------*........*--------*
        if (A2 > B1) or (A1 > B2):
            return None

        # Get intersection segment.
        # A1       A2       B1       B2
        # *--------*========*--------*
        # A2       A1       B2       B1
        # *--------*========*--------*
        A, B = max(A1, A2), min(B1, B2)

        # Intersection can be segment or single point.
        if A == B:
            return A
        else:
            Segment(A, B)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_plane(s, pl):
        """
        Intersection of segment and plane.

        Parameters
        ----------
        s : Segment
            Segment.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_triangle(s, t):
        """
        Intersection of segment and triangle.

        Parameters
        ----------
        s : Segment
            Segment.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_point(pl, p):
        """
        Intersection of plane and segment.

        Parameters
        ----------
        pl : Plane
            Plane.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on plane.
        """

        return Intersection.point_plane(p, pl)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_line(pl, ln):
        """
        Intersection of plane with line.

        Parameters
        ----------
        pl : Plane
            Plane.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Line
            Line is in plane.
        """

        return Intersection.line_plane(ln, pl)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_segment(pl, s):
        """
        Intersection of plane and segment.

        Parameters
        ----------
        pl : Plane
            Plane.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment is in plane.
        """

        return Intersection.segment_plane(s, pl)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_plane(pl1, pl2):
        """
        Intersection of two planes.

        Parameters
        ----------
        pl1 : Plane
            First plane.
        pl2 : Plane
            Second plane.

        Returns
        -------
        None
            Parallel planes.
        Line
            Intersect planes.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_triangle(pl, t):
        """
        Intersection of plane and triangle.

        Parameters
        ----------
        pl : Plane
            Plane.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment intersection.
        Triangle
            Triangle is in plane.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_point(t, p):
        """
        Intersection of triangle with point.

        Parameters
        ----------
        t : Triangle
            Triangle.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is in triangle.
        """

        return Intersection.point_triangle(p, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_line(t, ln):
        """
        Intersection of triangle and line.

        Parameters
        ----------
        t : Triangle
            Triangle.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Single point intersection.
        Segment
            Intersection by segment.
        """

        return Intersection.line_triangle(ln, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_segment(t, s):
        """
        Intersection of triangle with segment.

        Parameters
        ----------
        t : Triangle
            Triangle.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        return Intersection.segment_triangle(s, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_plane(t, pl):
        """
        Intersection triangle with plane.

        Parameters
        ----------
        t : Triangle
            Triangle.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment intersection.
        Triangle
            Triangle is in plane.
        """

        return Intersection.plane_triangle(pl, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_triangle(t1, t2):
        """
        Intersection of two triangles.

        Parameters
        ----------
        t1 : Triangle
            First triangle.
        t2 : Triangle
            Second triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        Triangle
            Triangle intersection.
        Points
            Convex figure of intersection
            4 - convex quadrangle,
            5 - convex pentagon,
            6 - convex hexagon.
        """

        raise Exception('not implemented')

#===================================================================================================

# Global objects.
O = Point(Fr(0), Fr(0), Fr(0))
X = Point(Fr(1), Fr(0), Fr(0))
Y = Point(Fr(0), Fr(1), Fr(0))
Z = Point(Fr(0), Fr(0), Fr(1))
OX = Line.from_points(O, X)
OY = Line.from_points(O, Y)
OZ = Line.from_points(O, Z)
XY = Line.from_points(X, Y)
OXY = Plane.from_points(O, X, Y)
OYZ = Plane.from_points(O, Y, Z)
OXZ = Plane.from_points(O, X, Z)
XYZ = Plane.from_points(X, Y, Z)

#===================================================================================================

def test():
    """
    Tests.
    """

    # Perpendicular planes.
    assert OXY.is_perpendicular_with_plane(OYZ)
    assert OXY.is_perpendicular_with_plane(OXZ)
    assert OYZ.is_perpendicular_with_plane(OXZ)

    # Intersect plane with line.
    assert XYZ.intersection_with_line(OX) == X
    assert XYZ.intersection_with_line(OY) == Y
    assert XYZ.intersection_with_line(OZ) == Z
    assert OXY.intersection_with_line(OX) == OX
    assert OYZ.intersection_with_line(OY) == OY
    assert OXZ.intersection_with_line(OZ) == OZ

    # Check equal lines.
    assert OX == Line.from_points(O, Point(Fr(2), Fr(0), Fr(0)))

    # Check point on segment.
    SOX = Segment(O, X)
    assert SOX.is_have_point(O)
    assert SOX.is_have_point(Point(Fr(1, 2), Fr(0), Fr(0)))
    assert not SOX.is_have_point(Y)
    assert not SOX.is_have_point(Point(Fr(3, 2), Fr(0), Fr(0)))

    # Check point on triangle.
    TXYZ = Triangle(X, Y, Z)
    assert TXYZ.is_have_point(Point(Fr(1, 3), Fr(1, 3), Fr(1, 3)))
    assert TXYZ.is_have_point(Point(Fr(2, 5), Fr(2, 5), Fr(1, 5)))
    assert not TXYZ.is_have_point(Point(Fr(2, 3), Fr(2, 3), Fr(-1, 3)))
    assert not TXYZ.is_have_point(Point(Fr(-1), Fr(-1), Fr(-1)))

    # Check lines intersections.
    assert OX.intersection_with_line(OY) == O
    assert OX.intersection_with_line(OZ) == O

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test()

#===================================================================================================
