import copy
import math
from fractions import Fraction as Fr
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
        x : number | Fraction
            x-coord.
        y : number | Fraction
            y-coord.
        """

        self.x = x
        self.y = y

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

        if size > 0:
            plt.scatter(self.x, self.y, color=color, s=size)

    #-----------------------------------------------------------------------------------------------

    def is_segment_end(self, s):
        """
        Check if point end of segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if point is end of segment, False - otherwise.
        """

        return (self == s.a) or (self == s.b)

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

        assert a != b
        self.a = a
        self.b = b
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
            Size of point.
        """

        x = [self.a.x, self.b.x]
        y = [self.a.y, self.b.y]
        if linewidth != '0':
            plt.plot(x, y, color=color, linewidth=linewidth)
        self.a.draw(plt, color=color, size=size)
        self.b.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def mod2(self):
        """
        Square of module.

        Returns
        -------
        number
            Square of length.
        """

        return (self.a.x - self.b.x)**2 + (self.a.y - self.b.y)**2

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
                # Two segments lie on one line.
                assert self.a < self.b
                assert s.a < s.b
                if self.b == s.a:
                    return self.b
                elif s.b == self.a:
                    return s.b
                elif (self.b < s.a) or (s.b < self.a):
                    return None
                elif self.b > s.a:
                    return Segment(s.a, self.b)
                elif s.b > self.a:
                    return Segment(self.a, s.b)
                else:
                    assert False
            else:
                return line.intersection_with_line(s_line)

    #-----------------------------------------------------------------------------------------------

    def sort_points(self):
        """
        Sort points.
        """

        if self.a > self.b:
            self.a, self.b = self.b, self.a

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
        Draw segments.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Line width.
        size : int
            Point size.
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

    def add_not_conflict(self, s):
        """
        Add new segment without conflicts.
        New segment can be added if:
        - it doesn't intersect any segment,
        - it has one common point with some segment.
        New segment can not be added if:
        - it intersect some segment and there is common point - inner point of one of them,
        - same segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if new segment added, False - otherwise.
        """

        for si in self.items:
            if si == s:
                return False
            else:
                intersec = si.intersection_with_segment(s)
                if isinstance(intersec, Point):
                    if (not intersec.is_segment_end(si)) and (not intersec.is_segment_end(s)):
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
                ss.append(Segment(ps[j], ps[j - 1]))

        # Set new set of segments.
        self.items = ss

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
            ps.add_unique(s.a)
            ps.add_unique(s.b)

        return ps

    #-----------------------------------------------------------------------------------------------

    def sort(self, fun):
        """
        Sort set of segments.

        Parameters
        ----------
        fun : fun
            Compare function.
        """

        self.items.sort(key=fun)

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
            return self
        else:
            a1, b1, c1 = self.a, self.b, self.c
            a2, b2, c2 = line.a, line.b, line.c
            if a1 == 0:
                #        b1 y + c1 = 0
                # a2 x + b2 y + c2 = 0
                y = -c1 / b1
                assert a2 != 0
                x = -(b2 * y + c2) / a2
                return Point(x, y)
            else:
                d = a1 * b2 - a2 * b1
                if d == 0:
                    assert False
                else:
                    y = (c1 * a2 - c2 * a1) / d
                    x = -(b1 * y + c1) / a1
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
            Points size.
        """

        self.ab.draw(plt, color=color, linewidth=linewidth)
        self.bc.draw(plt, color=color, linewidth=linewidth)
        self.ac.draw(plt, color=color, linewidth=linewidth)
        self.a.draw(plt, color=color, size=size)
        self.b.draw(plt, color=color, size=size)
        self.c.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def triangulation_segments(self, ps, ss):
        """
        Get triangulation segments.
        New segments must not intersect already existing segments.

        Parameters
        ----------
        ps : Points
            Given set of points.
        ss : Segments
            Given set of already existing segments.

        Returns
        -------
        Segments
            New set of segments.
        """

        # Create list of all possible segments.
        segments_all = Segments()
        pn = ps.count()
        for i in range(pn):
            for j in range(i + 1, pn):
                segments_all.add_unique(Segment(ps[i], ps[j]))
        segments_all.add_unique(self.ab)
        segments_all.add_unique(self.ac)
        segments_all.add_unique(self.bc)
        for i in range(pn):
            segments_all.add_unique(Segment(ps[i], self.a))
            segments_all.add_unique(Segment(ps[i], self.b))
            segments_all.add_unique(Segment(ps[i], self.c))
        segments_all.sort(fun=lambda s: s.mod2())

        # Create triangulation segments.
        segments_tri = Segments()
        for s in segments_all.items:
            segments_tri.add_not_conflict(s)

        return segments_tri

#===================================================================================================

def test():
    """
    Tests.
    """

    # Test lines.
    p00 = Point(Fr(0), Fr(0))
    p10 = Point(Fr(1), Fr(0))
    p01 = Point(Fr(0), Fr(1))
    line = Line.from_points(p00, p10)
    assert (line.a == 0) and (line.b == 1) and (line.c == 0)
    line = Line.from_points(p00, p01)
    assert (line.a == 1) and (line.b == 0) and (line.c == 0)
    line = Line.from_points(p10, p01)
    assert (line.a == 1) and (line.b == 1) and (line.c == -1)

    # Test segments intersections.
    s = Segment(Point(Fr(1), Fr(1)), Point(Fr(4), Fr(4)))
    s1 = Segment(Point(Fr(4), Fr(1)), Point(Fr(1), Fr(4)))
    s2 = Segment(Point(Fr(6), Fr(1)), Point(Fr(4), Fr(3)))
    assert s.is_intersects_with_segment(s1)
    assert s.intersection_with_segment(s1) == Point(Fr(5, 2), Fr(5, 2))
    assert not s.is_intersects_with_segment(s2)
    assert s.intersection_with_segment(s2) is None

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test()
    t = Triangle(Point(Fr(1), Fr(1)), Point(Fr(11), Fr(1)), Point(Fr(6), Fr(9)))
    ss = Segments()
    ss.add(Segment(Point(Fr(3), Fr(2)), Point(Fr(8), Fr(4))))
    ss.add(Segment(Point(Fr(6), Fr(2)), Point(Fr(6), Fr(7))))
    ss.add(Segment(Point(Fr(8), Fr(1)), Point(Fr(7), Fr(5))))
    ss.add(Segment(Point(Fr(4), Fr(4)), Point(Fr(9), Fr(3))))
    ss.add(Segment(Point(Fr(5), Fr(6)), Point(Fr(7), Fr(7))))
    ss.split_by_intersections()
    tri_ss = t.triangulation_segments(ss.points(), ss)
    tri_ss.draw(plt, color='blue', linewidth='1', size=0)
    t.draw(plt, color='black', linewidth='3', size=50)
    ss.draw(plt, color='black', linewidth='3', size=50)
    plt.show()

#===================================================================================================
