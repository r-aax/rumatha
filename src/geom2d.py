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
        x : any
            x-coord.
        y : any
            y-coord.
        """

        self.x = x
        self.y = y

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

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Line ({self.a}) x + ({self.b}) y + ({self.c})'

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

if __name__ == '__main__':
    # Initial points.
    p00 = Point(Fraction(0, 1), Fraction(0, 1))
    p10 = Point(Fraction(1, 1), Fraction(0, 1))
    p01 = Point(Fraction(0, 1), Fraction(1, 1))
    t = Triangle(p00, p10, p01)
    # Test lines.
    line = Line.from_points(p00, p10)
    assert (line.a == 0) and (line.b == 1) and (line.c == 0)
    line = Line.from_points(p00, p01)
    assert (line.a == 1) and (line.b == 0) and (line.c == 0)
    line = Line.from_points(p10, p01)
    assert (line.a == 1) and (line.b == 1) and (line.c == -1)
    #
    print(t)
    t.draw(plt)
    plt.show()
    pass

#===================================================================================================
