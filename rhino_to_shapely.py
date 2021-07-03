import rhino3dm as rh
from shapely.geometry import MultiLineString, asLineString, asPoint
from shapely.ops import polygonize
import numpy as np
import os.path

class CoordTransform:
    """Class that is responsible for the tranformations between a
    3d environment (x,y,z) and a 2d environment (x',y'). The 2d plane is derived from
    two vectors that define the 2d plane.

    Parameters
    ----------
    vec1 : ndarray
        A 3d vector (x,y,z) in the prefered 2d plane. This is normalised and used as e1.
    vec2 : ndarray
        Another 3d vector (x,y,z) in the prefered 2d plane. 

    Attributes
    ----------
    plane_normal : ndarray
        The normal of the plane in (x,y,z).

    Methods
    -------
    transform(pnts) :
        Transforms a coordinate from the (x,y,z) into (x',y')
    """
    def __init__(self, vec1, vec2):
        """Constructor

        Parameters
        ----------
        vec1 : ndarray
            A 3d vector (x,y,z) in the prefered 2d plane. This is normalised and used as e1.
        vec2 : ndarray
            Another 3d vector (x,y,z) in the prefered 2d plane. 
        """
        self._e1 = vec1 / np.linalg.norm(vec1)
        e3 = np.cross(vec1, vec2) 
        self._e3 = e3 / np.linalg.norm(e3)
        self._e2 = np.cross(self._e3,self._e1)
        self.Tinv = np.array([self._e1,self._e2,self._e3]).T
        self.T = np.linalg.inv(self.Tinv)
    
    def transform(self, pnts):
        """Transforms a coordinate from the (x,y,z) into (x',y')

        Parameters
        ----------
        pnts : ndarray
            Points in (x,y,z)

        Returns
        -------
        ndarray
            Points in (x',y')
        """
        pnts_prime = self.T.dot(pnts)
        if pnts.ndim==1: return pnts_prime[0:2]
        return pnts_prime[0:2,:]

    @property
    def plane_normal(self):
        return self._e3

class RhImporter:
    """Import geometric objects from rhino.

    Parameters
    ----------
    file_name : string
        The file name of a *.3dm file.
    brep : rhino3dm.Brep
        A rhino3dm brep

    Methods
    -------
    from_file(file_name) :
        Classmethod to create the object from a file.
    from_serialzed_brep(s_brep):
        Classmethod to create the object from a serialized brep object.
    get_planer_surface(refine_num, vec1, vec2, plane_distance, project, parallel)
        Generator that returns that single surface planer breps as shapely polygons.
    """
    def __init__(self, **kwarg):
        if "file_name" in kwarg:
            self._file_name = kwarg["file_name"]
            model = rh.File3dm.Read(self._file_name)
            self._process_objects(model.Objects)
        elif "brep" in kwarg:
            self._brep = [kwarg["brep"]]
    
    @classmethod
    def from_file(cls, file_name):
        """Class method to read from a Rhino file.

        Parameters
        ----------
        file_name : string
            File name of (or path to) a Rhino file.

        Returns
        -------
        RhImporter

        Raises
        ------
        ValueError
            If file extention is not ".3dm" or the file does not exit.
        """
        if not cls._validate_file_name(file_name):
            raise ValueError("File does not exist or is not a Rhino file.")
        return cls(file_name=file_name)

    @classmethod
    def from_serialzed_brep(cls, s_brep):
        """Class method to read from a serialized brep object.

        Parameters
        ----------
        s_brep : string
            Serialization of a brep

        Returns
        -------
        RhImporter

        Raises
        ------
        ValueError
            If the provided serialized object is not a single surface planer brep
        """
        brep = rh.CommonObject.Decode(s_brep)
        if not cls._validate_brep(brep): raise ValueError("Data is not a single surface planer brep.")
        return cls(brep=brep)

    @staticmethod
    def _validate_file_name(file_name):
        """Performs checks on a file name.

        Parameters
        ----------
        file_name : string
            The file name

        Returns
        -------
        Boolean
        """
        valid = False
        if os.path.isfile(file_name):
            valid = file_name.endswith('.3dm')
        return valid

    @staticmethod
    def _validate_brep(geom):
        """Performs checks on a brep.

        Parameters
        ----------
        geom : rhino3dm.GeometryBase
            Geometry to test for brep properties.

        Returns
        -------
        Boolean
        """
        valid = False
        if isinstance(geom,rh.Brep):
                if len(geom.Surfaces)==1:
                    valid = geom.Surfaces[0].IsPlanar()
        return valid

    @staticmethod
    def _validate_surface(geom):
        """Performs checks on a surface.

        Parameters
        ----------
        geom : rhino3dm.GeometryBase
            Geometry to test for surface properties.

        Returns
        -------
        Boolean
        """
        valid = False
        if isinstance(geom, rh.Surface):
            valid = geom.IsPlanar()
        return valid

    @staticmethod
    def _validate_curve(geom):
        """Performs checks on a curve.

        Parameters
        ----------
        geom : rhino3dm.GeometryBase
            Geometry to test for curve properties.

        Returns
        -------
        Boolean
        """
        return isinstance(geom, rh.Curve)

    @staticmethod
    def _validate_point(geom):
        """Performs checks on a point.

        Parameters
        ----------
        geom : rhino3dm.GeometryBase
            Geometry to test for point properties.

        Returns
        -------
        Boolean
        """
        return isinstance(geom, (rh.Point3d, rh.Point, rh.Point2d))

    def _process_objects(self, objects):
        """Process the file object table.

        Parameters
        ----------
        objects : rhino3dm.File3dmObjectTable
            The object table for the file.
        """
        self._brep = []
        self._surfaces = []
        self._curve = []
        self._point = []
        for obj in objects:
            if self._validate_brep(obj.Geometry): self._brep.append(obj.Geometry)
            elif self._validate_surface(obj.Geometry): self._surfaces.append(obj.Geometry)
            elif self._validate_curve(obj.Geometry): self._curve.append(obj.Geometry)
            elif self._validate_point(obj.Geometry): self._point.append(obj.Geometry)
            
    def get_planer_brep(self, refine_num, vec1=np.array([1,0,0]), vec2=np.array([0,1,0]), plane_distance=0.0, project=True, parallel=False):
        """Get all the single surface planer breps as Shapely polygons.
        Two vectors `vec1` and `vec2` describe the Shapely plane, with coordinates (x',y').
        The breps coordinates (x,y,z) are projected onto (x',y').
        Options to filter breps are provided:
        * only breps that are parallel to the Shapely plane
        * only breps in the Shapely plane
 
        Parameters
        ----------
        refine_num : integer
            Bézier curve interpolation number. In Rhino a surface's edges are nurb based curves.
            Shapely does not support nurbs, so the individual Bézier curves are interpolated using straight lines.
            This parameter sets the number of straight lines used in the interpolation.
        vec1 : numpy array, optional
            A 3d vector in the Shapely plane. Rhino is a 3D geometry environment.
            Shapely is a 2D geometric library.
            Thus a 2D plane needs to be defined in Rhino that represents the Shapely coordinate system.
            `vec1` represents the 1st vector of this plane. It will be used as Shapely's x direction.
        vec2 : numpy array, optional
            Continuing from `vec1`, `vec2` is another vector to define the Shapely plane.
            It must not be [0,0,0] and it's only requirement is that it is any vector in the Shapely plane (but not equal to `vec1`).
        plane_distance : float, optional
            The distance to the Shapely plane. 
        project : Boolean, optional
            Controls if the breps are projected onto the plane in the direction of the Shapley plane's normal.
        parallel : Boolean, optional
            Controls if only the rhino surfaces that have the same normal as the Shapely plane are yielded.
            If true, all non parallel surfaces are filtered out.
    
        Yields
        -------
        Shapely polygon.
            Shapely polygons with a coordinate system defined by vec1, and vec2.
    
        Raises
        ------
        ValueError
            If vec1 is not of ndim=1 and size=3
        ValueError
            If vec2 is not of ndim=1 and size=3
        ValueError
            If vec1 == vec2
        ValueError
            If vec1 or vec2 are the origin
        ValueError
            If both project and parallel are false. 
        """

        if not (vec1.ndim == 1 and vec1.size == 3):
            raise ValueError("vec1 is a numpy vector in 3d")
        if not (vec2.ndim == 1 and vec2.size == 3):
            raise ValueError("vec2 is a numpy vector in 3d")
        if not project and not parallel:
            raise ValueError("No surface meets this criteria, a surface that is not parallel and is not projected. This would just be the intersction of the plane and the surface (i.e. a line).")
        if (vec1==vec2).all():
            raise ValueError("vec2 must be different from vec1.")
        if (vec1==0).all() or (vec2==0).all():
            raise ValueError("The vectors must not be the origin.")

        ct = CoordTransform(vec1,vec2)

        def validation_factory():
            if project and parallel: return lambda normal, *args: (ct.plane_normal==np.array([normal.X,normal.Y,normal.Z])).all() 
            if project and not parallel: return lambda *args: True
            if not project and parallel: 
                return lambda normal, orig: (ct.plane_normal==np.array([normal.X,normal.Y,normal.Z])).all() and \
                    normal.X*orig.X + normal.Y*orig.Y + normal.Z*orig.Z == plane_distance
            
        validation = validation_factory()

        for surf in self._brep:
            _, frame = surf.Surfaces[0].FrameAt(0,0)
            if validation(frame.ZAxis, frame.Origin):
                rh_curvs = []
                for ii in range(0,len(surf.Edges)):
                    rh_curvs.append(RhCurv(surf.Edges[ii]))
                for rh_curv in rh_curvs:
                    if not rh_curv.is_line(): 
                        rh_curv.refine(refine_num)
                ml = MultiLineString([rc.get_shapely_line(ct.transform) for rc in rh_curvs])
                pgs = list(polygonize(ml))
                if len(pgs)>0: yield pgs[0]
    
    def get_curves(self, refine_num, vec1=np.array([1,0,0]), vec2=np.array([0,1,0]), plane_distance=0, project=True, parallel=False):
        if not (vec1.ndim == 1 and vec1.size == 3):
            raise ValueError("vec1 is a numpy vector in 3d")
        if not (vec2.ndim == 1 and vec2.size == 3):
            raise ValueError("vec2 is a numpy vector in 3d")
        if not project and not parallel:
            raise ValueError("No surface meets this criteria, a surface that is not parallel and is not projected. This would just be the intersction of the plane and the surface (i.e. a line).")
        if (vec1==vec2).all():
            raise ValueError("vec2 must be different from vec1.")
        if (vec1==0).all() or (vec2==0).all():
            raise ValueError("The vectors must not be the origin.")

        ct = CoordTransform(vec1,vec2)

        def validation_factory():
            if project and parallel: 
                # 1)check it is planer 2) check 2 points have the same distance value 
                return lambda planer, *args: planer and ct.plane_normal.dot(args[0]) - ct.plane_normal.dot(args[1]) == 0
            if project and not parallel: return lambda *args: True
            if not project and parallel: 
                # 1)check it is planer 2) check a point is in the plane 
                return lambda planer, *args: planer and ct.plane_normal.dot(args[0]) == plane_distance
            
        validation = validation_factory()

        for curve in self._curve:
            curve_w = RhCurv(curve)
            if validation(curve_w.is_planer, *curve_w.get_greville_points):
                if not curve_w.is_line(): 
                    curve_w.refine(refine_num)
                ls = curve_w.get_shapely_line(ct.transform)
                yield ls

    def get_points(self, vec1=np.array([1,0,0]), vec2=np.array([0,1,0]), plane_distance=0, project=True):
        """Get all the rhino points as Shapely points.
        Two vectors `vec1` and `vec2` describe the Shapely plane, with coordinates (x',y').
        The point coordinates (x,y,z) are projected onto (x',y').
        Options to filter breps are provided:
        * only points in the Shapely plane
 
        Parameters
        ----------
        vec1 : numpy array, optional
            A 3d vector in the Shapely plane. Rhino is a 3D geometry environment.
            Shapely is a 2D geometric library.
            Thus a 2D plane needs to be defined in Rhino that represents the Shapely coordinate system.
            `vec1` represents the 1st vector of this plane. It will be used as Shapely's x direction.
        vec2 : numpy array, optional
            Continuing from `vec1`, `vec2` is another vector to define the Shapely plane.
            It must not be [0,0,0] and it's only requirement is that it is any vector in the Shapely plane (but not equal to `vec1`).
        plane_distance : float, optional
            The distance to the Shapely plane.
        project : Boolean, optional
            Controls if the points are projected onto the plane in the direction of the Shapley plane's normal.
    
        Yields
        -------
        Shapely point.
            Shapely points with a coordinate system defined by vec1, and vec2.
    
        Raises
        ------
        ValueError
            If vec1 is not of ndim=1 and size=3
        ValueError
            If vec2 is not of ndim=1 and size=3
        ValueError
            If vec1 == vec2
        ValueError
            If vec1 or vec2 are the origin
        """

        if not (vec1.ndim == 1 and vec1.size == 3):
            raise ValueError("vec1 is a numpy vector in 3d")
        if not (vec2.ndim == 1 and vec2.size == 3):
            raise ValueError("vec2 is a numpy vector in 3d")
        if (vec1==vec2).all():
            raise ValueError("vec2 must be different from vec1.")
        if (vec1==0).all() or (vec2==0).all():
            raise ValueError("The vectors must not be the origin.")

        ct = CoordTransform(vec1,vec2)

        def validation_factory():
            if project: return lambda *args: True
            if not project: 
                return lambda pnt: ct.plane_normal.dot(pnt) == -plane_distance
        
        validation = validation_factory()
        for pnt in self._point:
            pnt_w = RhPnt(pnt)
            if validation(pnt_w.as_numpy):
                yield pnt_w.get_shapely_point(ct.transform)
        
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