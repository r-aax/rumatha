"""
3D geometry in rational coordinates.
"""

from fractions import Fraction as Fr
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

        line = Line.from_segment(self)

        return line.is_have_point(p) and p.is_between(self.A, self.B)

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

        assert (m != 0) or (n != 0) or (p != 0)

        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.m = m
        self.n = n
        self.p = p

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

    @staticmethod
    def from_segment(s):
        """
        Constructor from segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        Line
            Constructed line.
        """

        return Line.from_points(s.A, s.B)

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
        x0, y0, z0, m, n, p = self.x0, self.y0, self.z0, self.m, self.n, self.p

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
        x0, y0, z0, m, n, p = line.x0, line.y0, line.z0, line.m, line.n, line.p

        # ax + by + cz + d = 0
        #   x = x0 + tm
        #   y = y0 + tn
        #   z = z0 + tp

        # One of m, n, p is not zero.

        if line.m != 0:

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

        elif line.n != 0:

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

            assert line.p != 0

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

#===================================================================================================

# Global objects.
O = Point(Fr(0), Fr(0), Fr(0))
X = Point(Fr(1), Fr(0), Fr(0))
Y = Point(Fr(0), Fr(1), Fr(0))
Z = Point(Fr(0), Fr(0), Fr(1))
OX = Line.from_points(O, X)
OY = Line.from_points(O, Y)
OZ = Line.from_points(O, Z)
OXY = Plane.from_points(O, X, Y)
OYZ = Plane.from_points(O, Y, Z)
OXZ = Plane.from_points(O, X, Z)
XYZ = Plane.from_points(X, Y, Z)

#===================================================================================================

def test():
    """
    Tests.
    """

    # Intersect plane with line.
    assert XYZ.intersection_with_line(OX) == X
    assert XYZ.intersection_with_line(OY) == Y
    assert XYZ.intersection_with_line(OZ) == Z
    assert OXY.intersection_with_line(OX) == OX
    assert OYZ.intersection_with_line(OY) == OY
    assert OXZ.intersection_with_line(OZ) == OZ

    # Check point on segment.
    SOX = Segment(O, X)
    assert SOX.is_have_point(O)
    assert SOX.is_have_point(Point(Fr(1, 2), Fr(0), Fr(0)))
    assert not SOX.is_have_point(Y)
    assert not SOX.is_have_point(Point(Fr(3, 2), Fr(0), Fr(0)))

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test()

#===================================================================================================
