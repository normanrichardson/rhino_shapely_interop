import rhino3dm as rh
from shapely.geometry import MultiLineString, asLineString
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
    transform_rh_sh(pnts) :
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
        return self.T.dot(pnt)

    @property
    def plane_normal(self):
        return self._e3

class RhImporter:
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
    rh.File3dm.Objects
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
        valid = False
        if isinstance(geom, rh.Curve):
            if valid:=geom.IsLinear(): 
                valid = True
            else:
                valid = geom.IsPlanar()
        return valid

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
        return isinstance(geom, rh.Point3d)

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
            
    def get_planer_surface(self, refine_num, vec1=np.array([1,0,0]), vec2=np.array([0,1,0]), plane_distance=0.0, project=True, parallel=False):
        """Get all the single surface planer breps.
       Two vectors `vec1` and `vec2` describe the shapely plane, with coordinates (x',y').
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
           Continuing from `vec1`, `vec2` is another vector to befine the Shapely plane.
           It must not be [0,0,0] and it's only requirement is that it is any vector in the Shapely plane (but not equal to `vec1`).
       plane_distance : float, optional
           The distance to the Shapely plane. If it is not provided, all geometry is projected onto the provided plane(via `vec1` and `vec2`).
           If it is provided, then only surfaces in the unique plane (defined by `vec1`, `vec2`, and `plane_distance`) are yielded.
       project : Boolean, optional
           Controls if the shapes are projected onto the plane in the direction of the shapley plane's normal.
       parallel : Boolean, optional
           Controls if only the rhino surface have the same normal as the Shapely plane are yielded.
           If true, all non parallel surfaces are filtered out.
 
       Yields
       -------
       Shapely Surface.
           Shapely surface with a coordinate system defined by the rhino surface's normal vector.
 
       Raises
       ------
       ValueError
           If a plane distance is provided, but a plane normal is not.
       """

        if not (vec1.ndim == 1 and vec1.size == 3):
            raise ValueError("vec1 is a numpy vector in 3d")
        if not (vec2.ndim == 1 and vec2.size == 3):
            raise ValueError("vec2 is a numpy vector in 3d")
        if not project and not parallel:
            raise ValueError("No surface meets this criteria, a surface that is not parallel and is not projected. This would just be the intersction of the plane and the surface (i.e. a line).")

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
    
    def get_planer_curves(self, vec1=np.array([1,0,0]), vec2=np.array([0,1,0]), plane_distance=0, project=True, parallel=False):
        pass

    def get_points(self, vec1=np.array([1,0,0]), vec2=np.array([0,1,0]), plane_distance=0, project=True):
        pass

class RhCurv:
    def __init__(self, curv):
        self._curv = curv
        self._nurb = curv.ToNurbsCurve()
        self._greville_points_param = [self._nurb.GrevilleParameter(idx) for idx in range(len(self._nurb.Points))]
        self._degree = self._nurb.Order-1
        self._greville_points_param_modif = self._greville_points_param
    def refine(self,num):
        gen_interv = ((self._greville_points_param[ii], self._greville_points_param[ii+1]) for ii in range(len(self._greville_points_param)-1))
        self._greville_points_param_modif = [self._greville_points_param[0]]
        for ii, jj in gen_interv:
            self._greville_points_param_modif += list(np.linspace(ii,jj,num+2)[1:])
    def get_shapely_line(self, transform):
        pnts = []
        for t in self._greville_points_param_modif:
            pnt = self._curv.PointAt(t)
            pnts.append([pnt.X, pnt.Y, pnt.Z])
        pnts_np = transform(np.array(pnts).T)[0:2,:].round(decimals=12)
        return asLineString(pnts_np.T)
    def is_line(self):
        return self._curv.IsLinear()
