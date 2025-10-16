"""
3D geometry with real coordinated.
"""

import geom1d

#===================================================================================================

class Point:
    """
    Point.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, x=0.0, y=0.0, z=0.0):
        """
        Constructor.

        Parameters
        ----------
        x : float
            X coordinate.
        y : float
            Y coordinate.
        z : float
            Z coordinate.
        """

        self.x = x
        self.y = y
        self.z = z

#===================================================================================================

class Triangle:
    """
    Triangle.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, A, B, C):
        """
        Constructor from points.

        Parameters
        ----------
        A : Point
            A point.
        B : Point
            B point.
        C : Point
            C Point.
        """

        self.A = A
        self.B = B
        self.C = C

#===================================================================================================

class Box:
    """
    Box.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, xlo, xhi, ylo, yhi, zlo, zhi):
        """
        Constructor.

        Parameters
        ----------
        xlo : float
            X lo value.
        xhi : float
            X hi value.
        ylo : float
            Y lo value.
        yhi : float
            Y hi value.
        zlo : float
            Z lo value.
        zhi : float
            Z hi value.
        """

        self.sx = geom1d.Segment(xlo, xhi)
        self.sy = geom1d.Segment(ylo, yhi)
        self.sz = geom1d.Segment(zlo, zhi)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_triangle(t):
        """
        Construct box from triangle.

        Parameters
        ----------
        t : Triangle
            Triangle.

        Returns
        -------
        Box
            Constructed box.
        """

        return Box(min([t.A.x, t.B.x, t.C.x]), max([t.A.x, t.B.x, t.C.x]),
                   min([t.A.y, t.B.y, t.C.y]), max([t.A.y, t.B.y, t.C.y]),
                   min([t.A.z, t.B.z, t.C.z]), max([t.A.z, t.B.z, t.C.z]))

#===================================================================================================

if __name__ == '__main__':
    pass

#===================================================================================================
