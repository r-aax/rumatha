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

    def __init__(self, x=0, y=0, z=0):
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

        return (self.x == p.x) and (self.y == p.y) and (self.z == z)

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

#===================================================================================================

def test():
    """
    Tests.
    """

    assert True

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test()

#===================================================================================================
