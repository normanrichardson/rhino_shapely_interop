import rhino3dm as rh
from shapely.geometry import MultiLineString, LineString
from shapely.ops import polygonize
import matplotlib.pyplot as plt
import numpy as np
import json

# rh.Surface.IsPlanar()

class RhFile:
    def __init__(self, file_name):
        self.validate(file_name)
        self._file_name = file_name
        model = rh.File3dm.Read(file_name)
        self._objects = model.Objects
        #print(self._objects[1].Geometry)
        self.process_objects()
    @staticmethod
    def validate(file_name):
        #does the file exist
        pass
    @staticmethod
    def validate_geom(self, geom):
        valid = False
        if isinstance(geom, rh.Surface):
            valid = geom.IsPlanar()
        # if isinstance(geom, rh.Curve):
        #     valid = geom.IsPlanar()
        #     valid = geom.IsClosed
        return valid
        # Is surface planer, ok get boundary
        # not done yet... What transformation gets it onto the x,y plane.
        # test plane
    def process_surface(self):
        pass
    def process_objects(self):
        self._surfaces = []
        # curves = {}
        for obj in self._objects:
            # geom = obj.Geometry
            # if isinstance(geom, rh.Curve):
            #     frame = geom.FrameAt(0)
            #     p_normal = frame[1].ZAxis
            #     if p_normal in curves: curves[p_normal].append(geom)
            #     else: curves[p_normal] = [geom]
            # else:
            self._surfaces.append(obj.Geometry)

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