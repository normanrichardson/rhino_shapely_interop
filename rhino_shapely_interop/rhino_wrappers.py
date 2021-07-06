from shapely.geometry import asLineString, asPoint
import numpy as np

class RhCurv:
    """Wrapper for a rhino curve.

    Parameters
    ----------
    curv :
        A rhino3dm.Curve

    Methods
    -------
    refine(num) :
        Refine the individual Bézier curves of the rhino.Curve
    get_shapely_line(transform) :
        Get the shapely line string for the rhino curve.
    is_line() :
        Is the rhino line a straight line
    """
    def __init__(self, curv):
        """Constructor

        Parameters
        ----------
        curv : rhino3dm.Curve
            A rhino3dm curve
        """
        self._curv = curv
        self._nurb = curv.ToNurbsCurve()
        self._greville_points_param = [self._nurb.GrevilleParameter(idx) for idx in range(len(self._nurb.Points))]
        self._degree = self._nurb.Order-1
        self._greville_points_param_modif = self._greville_points_param

    def refine(self,num):
        """Refine the individual Bézier curves of the rhino.Curve

        Parameters
        ----------
        num : integer
            Number of refinements
        """
        gen_interv = ((self._greville_points_param[ii], self._greville_points_param[ii+1]) for ii in range(len(self._greville_points_param)-1))
        self._greville_points_param_modif = [self._greville_points_param[0]]
        for ii, jj in gen_interv:
            self._greville_points_param_modif += list(np.linspace(ii,jj,num+2)[1:])

    def get_shapely_line(self, transform):
        """Get the shapely line string for the rhino curve.

        Parameters
        ----------
        transform : func
            A function that transforms (3,n) ndarray into a new coordinate system.

        Returns
        -------
        Shapely.Geometry.LineString
            The discretized shapely representation of a rhino curve.
        """
        pnts = []
        for t in self._greville_points_param_modif:
            pnt = self._curv.PointAt(t)
            pnts.append([pnt.X, pnt.Y, pnt.Z])
        pnts_np = transform(np.array(pnts).T).round(decimals=12)
        return asLineString(pnts_np.T)
    
    def is_line(self):
        """Is the rhino line a straight line

        Returns
        -------
        Boolean
        """
        return self._curv.IsLinear()

    @property
    def get_greville_points(self):
        pnts = []
        for t in self._greville_points_param:
            pnt = self._curv.PointAt(t)
            pnts.append(np.array([pnt.X, pnt.Y, pnt.Z]))
        return pnts
    
    @property
    def is_planer(self):
        return self._curv.IsPlanar()

class RhPnt:
    def __init__(self, pnt):
        """Wrapper for a rhino point.

        Parameters
        ----------
        pnt : RhinoPoint3d or RhinoPoint2d
        """
        self._pnt = pnt
        try:
            self._pnt_np = np.array([pnt.X, pnt.Y, pnt.Z])
        except:
            try:
                self._pnt_np = np.array([pnt.X, pnt.Y, 0])
            except:
                self._pnt_np = np.array([pnt.Location.X, pnt.Location.Y, pnt.Location.Z])

    def get_shapely_point(self, transform):
        """Get the shapely point string for the rhino point.

        Parameters
        ----------
        transform : func
            A function that transforms (3,n) ndarray into a new coordinate system.

        Returns
        -------
        Shapely.Geometry.PointString
            The shapely representation of a rhino point.
        """
        pnts_np = transform(np.array(self._pnt_np).T).round(decimals=12)
        return asPoint(pnts_np)

    @property
    def as_numpy(self):
        return self._pnt_np