import copy
import math
from fractions import Fraction
import matplotlib.pyplot as plt

#===================================================================================================

class Point:
    """
    Point with two coordinates.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, x, y):
        """
        Constructor.

        Parameters
        ----------
        x : number
            x-coord.
        y : number
            y-coord.
        """

        self.x = x
        self.y = y

    #-----------------------------------------------------------------------------------------------

    def __copy__(self):
        """
        Copy constructor.

        Returns
        -------
        Point
            Copy point.
        """

        return Point(self.x, self.y)

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, p):
        """
        Check equal to another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if equal to another point, False - otherwise.
        """

        return (self.x == p.x) and (self.y == p.y)

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
            True - it not equal to another point, False - otherwise.
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
            True - if self >= p, False - otherwise.
        """

        if self.x > p.x:
            return True
        elif self.x < p.x:
            return False
        else:
            return self.y >= p.y

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
            True - if self > p, False - otherwise.
        """

        if self.x > p.x:
            return True
        elif self.x < p.x:
            return False
        else:
            return self.y > p.y

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
            True - if self <= p, False - otherwise.
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
            True - if self < p, False - otherwise.
        """

        return not (self >= p)

    #-----------------------------------------------------------------------------------------------

    def  __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'P({self.x}, {self.y})'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        """

        plt.scatter(self.x, self.y, color='black', s=20)

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

    def draw(self, plt):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        """

        for p in self.items:
            p.draw(plt)

    #-----------------------------------------------------------------------------------------------

    def count(self):
        """
        Count of items.

        Returns
        -------
        int
            Count of items.
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
        Add unique point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point was added, False - otherwise.
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
    Segment.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, a, b):
        """
        Constructor.

        Parameters
        ----------
        a : Point
            A point.
        b : Point
            B point.
        """

        self.a = a
        self.b = b

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
            True - if equal to another segment, False - otherwise.
        """

        return ((self.a == s.a) and (self.b == s.b)) or ((self.a == s.b) and (self.b == s.a))

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
            True - if not equal to another segment, False - otherwise.
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

        return f'Segm({self.a}, {self.b})'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        """

        x = [self.a.x, self.b.x]
        y = [self.a.y, self.b.y]
        plt.plot(x, y, color='black', linewidth='2')
        self.a.draw(plt)
        self.b.draw(plt)

    #-----------------------------------------------------------------------------------------------

    def is_intersects_with_line(self, line):
        """
        Check if segment intersects with line.
        There is intersection if segment ends lie in the different sides from the line.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if intersects, False - otherwise.
        """

        return line.is_intersects_with_segment(self)

    #-----------------------------------------------------------------------------------------------

    def is_intersects_with_segment(self, s):
        """
        Check if segment intersects with another segment.

        Parameters
        ----------
        s : Segment.
            Segment.

        Returns
        -------
        bool
            True - if intersects, False - otherwise.
        """

        line = Line.from_points(self.a, self.b)
        s_line = Line.from_points(s.a, s.b)

        return self.is_intersects_with_line(s_line) and s.is_intersects_with_line(line)

    #-----------------------------------------------------------------------------------------------

    def intersection_with_segment(self, s):
        """
        Find intersection with another segment.

        Parameters
        ----------
        s : Segment.

        Returns
        -------
        None
            If there is no intersetion.
        Point
            If there is intersection on one point.
        Segment
            If there is intersection on segment.
        """

        if not self.is_intersects_with_segment(s):
            return None
        else:
            line = Line.from_segment(self)
            s_line = Line.from_segment(s)
            if line == s_line:
                # Two segment intersect on segment.
                assert False
            else:
                return line.intersection_with_line(s_line)

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

    def draw(self, plt):
        """
        Draw segments.

        Parameters
        ----------
        plt : Plot
            Plot.
        """

        for s in self.items:
            s.draw(plt)

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
        Add new unique segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if new segment was added, False - otherwise.
        """

        for si in self.items:
            if si == s:
                return False

        self.add(s)

        return True

    #-----------------------------------------------------------------------------------------------

    def split_by_intersections(self):
        """
        Find all intersection points and split all segments.
        """

        # Allocate data for points.
        n = self.count()
        pss = []
        for _ in range(n):
            pss.append(Points())

        # Find intersections of segments pairs.
        for i in range(n):
            for j in range(i + 1, n):
                si, sj = self.items[i], self.items[j]
                intersec = si.intersection_with_segment(sj)
                if isinstance(intersec, Point):
                    pss[i].add_unique(intersec)
                    pss[j].add_unique(intersec)
                else:
                    # Check that we have not intersection by segment.
                    assert intersec is None

        # Add segment ends to points sets and sort all points sets.
        for i in range(n):
            s = self.items[i]
            ps = pss[i]
            ps.add_unique(s.a)
            ps.add_unique(s.b)
            ps.sort()

        # Create new list of segments.
        ss = []
        for i in range(n):
            ps = pss[i]
            k = ps.count()
            for j in range(1, k):
                ss.append(Segment(copy.copy(ps[j]), copy.copy(ps[j - 1])))

        # Set new set of segments.
        self.items = ss

#===================================================================================================

class Line:
    """
    Line.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, a, b, c):
        """
        Constructor from coefficients.

        Parameters
        ----------
        a : any
            a coefficient.
        b : any
            b coefficient.
        c : any
            c coefficient.
        """

        self.a = a
        self.b = b
        self.c = c

    #-----------------------------------------------------------------------------------------------

    def __copy__(self):
        """
        Line copy.

        Returns
        -------
        Line
            New line.
        """

        return Line(self.a, self.b, self.c)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_points(a, b):
        """
        Create line from two points.

        Parameters
        ----------
        a : Point
            A point.
        b : Point
            B point.

        Returns
        -------
        Line
            Result line.
        """

        x1, y1, x2, y2 = a.x, a.y, b.x, b.y

        if y1 == y2:
            if x1 == x2:
                raise Exception(f'Line from a = {a}, b = {b} can not be constructed')
            else:
                # coefficients must be initialized with int values
                # (for compatibility with both floats and Fractions)
                return Line(a=0, b=1, c=-y1)
        else:
            ka = 1
            kb = (x2 - x1) / (y1 - y2)
            kc = (x1 * y2 - x2 * y1) / (y1 - y2)
            return Line(a=ka, b=kb, c=kc)

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
            Result line.
        """

        return Line.from_points(s.a, s.b)

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, line):
        """
        Check equal to another line.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if equal to another line, False - otherwise.
        """

        return (self.a == line.a) and (self.b == line.b) and (self.c == line.c)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, line):
        """
        Check not equal to another line.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if not equal to another line, False - otherwise.
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

        return f'Line ({self.a}) x + ({self.b}) y + ({self.c})'

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
        number
            Value of the point according to the line.
        """

        return self.a * p.x + self.b * p.y + self.c

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if line has  point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if line has point, False - otherwise.
        """

        return self.val(p) == 0

    #-----------------------------------------------------------------------------------------------

    def is_two_points_strong_on_one_side(self, p1, p2):
        """
        Check if two points strong on one side.
        These two points don't lie on the line.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point
            Second point.

        Returns
        -------
        bool
            True - if two points lie on one side of the line, False - otherwise.
        """

        return self.val(p1) * self.val(p2) > 0

    #-----------------------------------------------------------------------------------------------

    def is_intersects_with_segment(self, s):
        """
        Check if line intersects with segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if there is intersection, False - otherwise.
        """

        return not self.is_two_points_strong_on_one_side(s.a, s.b)

    #-----------------------------------------------------------------------------------------------

    def is_intersects_with_line(self, line):
        """
        Check if intersection with another line.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if there is intersection, False - otherwise.
        """

        if (self.a == line.a) and (self.b == line.b):
            return self.c == line.c
        else:
            return True

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
            If there is one point of intersection.
        Line
            If it is the same line.
        """

        if self == line:
            # Create new line.
            return copy.copy(self)
        else:
            a1, b1, c1 = self.a, self.b, self.c
            a2, b2, c2 = line.a, line.b, line.c
            if a1 == 0:
                assert False
            else:
                d = a1 * b2 - a2 * b1
                if d == 0:
                    assert False
                else:
                    y = (c1 * a2 - c2 * a1) / d
                    x = (-(b1 * y + c1)) / a1
                    return Point(x, y)

#===================================================================================================

class Triangle:
    """
    Triangle.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, a, b, c):
        """
        Constructor.

        Parameters
        ----------
        a : Point
            A point.
        b : Point
            B point.
        c : Point
            C point.
        """

        self.a = a
        self.b = b
        self.c = c

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Tri({self.a}, {self.b}, {self.c})'

    #-----------------------------------------------------------------------------------------------

    @property
    def ab(self):
        """
        AB side.

        Returns
        -------
        Segment
            AB side.
        """

        return Segment(self.a, self.b)

    #-----------------------------------------------------------------------------------------------

    @property
    def bc(self):
        """
        BC side.

        Returns
        -------
        Segment
            BC side.
        """

        return Segment(self.b, self.c)

    #-----------------------------------------------------------------------------------------------

    @property
    def ac(self):
        """
        AC side.

        Returns
        -------
        Segment
            AC side.
        """

        return Segment(self.a, self.c)

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        """

        self.a.draw(plt)
        self.b.draw(plt)
        self.c.draw(plt)
        self.ab.draw(plt)
        self.bc.draw(plt)
        self.ac.draw(plt)

#===================================================================================================

def test():
    """
    Tests.
    """

    # Test lines.
    p00 = Point(Fraction(0, 1), Fraction(0, 1))
    p10 = Point(Fraction(1, 1), Fraction(0, 1))
    p01 = Point(Fraction(0, 1), Fraction(1, 1))
    line = Line.from_points(p00, p10)
    assert (line.a == 0) and (line.b == 1) and (line.c == 0)
    line = Line.from_points(p00, p01)
    assert (line.a == 1) and (line.b == 0) and (line.c == 0)
    line = Line.from_points(p10, p01)
    assert (line.a == 1) and (line.b == 1) and (line.c == -1)

    # Test segments intersections.
    s = Segment(Point(1, 1), Point(4, 4))
    s1 = Segment(Point(4, 1), Point(1, 4))
    s2 = Segment(Point(6, 1), Point(4, 3))
    assert s.is_intersects_with_segment(s1)
    assert s.intersection_with_segment(s1) == Point(Fraction(5, 2), Fraction(5, 2))
    assert not s.is_intersects_with_segment(s2)
    assert s.intersection_with_segment(s2) is None

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test()
    t = Triangle(Point(1, 1), Point(11, 1), Point(6, 9))
    ss = Segments()
    ss.add(Segment(Point(3, 2), Point(8, 4)))
    ss.add(Segment(Point(6, 2), Point(6, 7)))
    ss.add(Segment(Point(8, 2), Point(7, 5)))
    ss.add(Segment(Point(4, 4), Point(9, 3)))
    ss.add(Segment(Point(5, 6), Point(7, 7)))
    ss.split_by_intersections()
    t.draw(plt)
    ss.draw(plt)
    print(ss.items)
    plt.show()

#===================================================================================================
