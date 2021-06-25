import rhino3dm as rh
from shapely.geometry import MultiLineString, LineString
from shapely.ops import polygonize
import matplotlib.pyplot as plt
import numpy as np
import json

class CoordTransform:
    def __init__(self, vec1, vec2) -> None:
        self._e1 = vec1 / np.linalg.norm(vec1)
        e3 = np.cross(vec1, vec2) 
        self._e3 = e3 / np.linalg.norm(e3)
        self._e2 = np.cross(self._e3,self._e1)
        self.Tinv = np.matrix([self._e1,self._e2,self._e3]).T
        self.T = np.linalg.inv(self.Tinv)
    
    @staticmethod
    def transform_rh_sh(self, pnt):
        return self.T.dot(pnt)

class RhImporter:
    def __init__(self, **kwarg):
        if "file_name" in kwarg:
            self._file_name = kwarg["file_name"]
            model = rh.File3dm.Read(self._file_name)
            self._process_objects(model.Objects)
        elif "surface" in kwarg:
            self._surfaces = kwarg["surface"]
    
    @classmethod
    def from_file(cls, file_name):
        if not cls._validate_file_name(file_name):
            raise ValueError("File name not valid.")
        return cls.__init__(file_name=file_name)

    @classmethod
    def from_serialzed_surface(cls, s_surface):
        # need some validation here...
        surface = rh.CommonObject.Decode(s_surface)
        if not cls._validate_surface(surface): raise ValueError("Data is not surface or the surface it is not planer")
        return cls.__init__(surface=surface)

    @staticmethod
    def _validate_file_name(file_name):
        #does the file exist
        pass

    @staticmethod
    def _validate_surface(self, geom):
        valid = False
        if isinstance(geom, rh.Surface):
            valid = geom.IsPlanar()
        return valid

    @staticmethod
    def _validate_curve(self, geom):
        valid = False
        if isinstance(geom, rh.Curve):
            if valid:=geom.IsLinear(): 
                valid = True
            else:
                valid = geom.IsPlanar()
        return valid

    @staticmethod
    def _validate_point(self, geom):
        return isinstance(geom, rh.Point3d)

    def _process_objects(self, objects):
        self._surfaces = []
        self._curve = []
        self._point = []
        for obj in objects:
            if self._validate_surface(obj): self._surfaces.append(obj)
            if self._validate_curve(obj): self._curve.append(obj)
            if self._validate_point(obj): self._point.append(obj)
            
    def get_planer_surface(self, refine_num, vec1=None, vec2=None, plane_distance=0, project=False):

        """Get all the planer surfaces. If a plane normal is provided then only surfaces with this normal are returned (nonunique plane).
        If both a plane normal and a plane distance is provied, then only surfaces in this unique plane are returned.

        Parameters
        ----------
        refine_num : integer
            Bézier curve interpolation number. In Rhino a surfaces edges are nurb based curves. 
            Shapely does not support nurbs, so the individule Bézier curves are interplated using straight lines.
            This parameter sets the number of staight lines used in the interpolation.
        vec1 : numpy array, optional
            A 3d vector in the Shapely plane. Rhino is a 3D geometry envionoment.
            Shapely is a 2D geometric library.
            Thus a 2D plane needs to be defined in Rhino that represents the shapely coordinate system.
            `vec1` represents the 1st vector of this plane. It will be used as shapely's x direction.
        vec2 : numpy array, optional
            Continuing from `vec1`, `vec2` is another vector to befine the shapely plane.
            It must not be [0,0,0] and its only requirments is that it is any vector in the shapely plane (but not equal to `vec1`).
        plane_distance : float, optional
            The distance to the shapely plane. If it is not provided, all geometry is projected onto the provided plane(via `vec1` and `vec2`).
            If it is provided, then only planar.

        Yields
        -------
        Shapely Surface.
            Shapely surface with coordinates system defined by the rhino surfaces normal vector.

        Raises
        ------
        ValueError
            If a plane distance is provided, but a plane normal is not.
        """
        
        # construct the transformation matrix

        
        #if (vec2 is not None and vec1 is None) or (vec1 is not None and vec2 is None): 
        #    raise ValueError("If the Shapely plane is not fully defined. Provide both vec1 and vec2")

        #def validation_factory():
            #if vec1 is None and plane_distance is None: return lambda *args: True
            #if vec1 is not None and plane_distance is None: return lambda normal, *args: plane_normal == normal
            #if plane_normal is None and plane_distance is None: 
            #    return lambda normal, orig: plane_normal == normal and \
            #        normal.X*orig.X + normal.Y*orig.Y + normal.Z*orig.Z == plane_distance
        validation = validation_factory()

        for surf in self._surfaces:
            frame = surf.FrameAt(0,0)
            if validation(frame.ZAxis, frame.Origin):
                rh_curvs = []
                for ii in range(0,len(surf.Edges)):
                    rh_curvs.append(RhCurv(surf.Edges[ii]))
                for rh_curv in rh_curvs:
                    if not rh_curv.is_line(): 
                        rh_curv.refine(refine_num)
                ml = MultiLineString([rc.get_shapely_line() for rc in rh_curvs])
                pgs = list(polygonize(ml))
                yield pgs[0]
    
    def get_all_planer_curves(self, plane_normal, plane_distance=None):
        pass

    def get_all_points(self, plane_normal, plane_distance):
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
    def get_shapely_line(self):
        line = []
        for t in self._greville_points_param_modif:
            pnt = self._curv.PointAt(t)
            line.append((pnt.X, pnt.Y, pnt.Z))
        return LineString(line)
    def is_line(self):
        return self._curv.IsLinear()

def convert_and_refine(rh_geom, refine_num):
    rh_curvs = []
    for ii in range(0,len(rh_geom.Edges)):
        rh_curvs.append(RhCurv(rh_geom.Edges[ii]))
    for rh_curv in rh_curvs:
        if not rh_curv.is_line(): 
            rh_curv.refine(refine_num)
    ml = MultiLineString([rc.get_shapely_line() for rc in rh_curvs])
    pgs = list(polygonize(ml))
    return pgs[0], pgs[1:]

refine = 0

with open('data',"r") as f:
    json_data = f.read()
brep_dict = json.loads(json_data)
brep = rh.CommonObject.Decode(brep_dict[2])

ploy1, holes1 = convert_and_refine(brep, refine)

[plt.plot(*ploy1.exterior.xy)]
[plt.plot(*ii.xy) for ii in ploy1.interiors]
plt.show()

RhFile('AAA.3dm')