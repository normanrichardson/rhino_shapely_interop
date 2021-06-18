import rhino3dm as rh
from shapely.geometry import MultiLineString
from shapely.ops import polygonize
import matplotlib.pyplot as plt
import numpy as np
import json

with open('data',"r") as f:
    brep_dict = json.load(f)

lines = []
brep = rh.CommonObject.Decode(brep_dict)
for ii in range(0,len(brep.Edges)):
    if brep.Edges[ii].IsLinear():
        s_p = brep.Edges[ii].PointAtStart
        e_p = brep.Edges[ii].PointAtEnd
        l_points = [(s_p.X, s_p.Y, s_p.Z), (e_p.X, e_p.Y, e_p.Z)]
        lines.append(l_points)
    else: 
        nurb = brep.Edges[ii].ToNurbsCurve()
        l_points = []
        for idx, point in enumerate(nurb.Points):
            if point.W==1: l_points.append((point.X, point.Y, point.Z))
        lines.append(l_points)

rh.NurbsCurve.MakePiecewiseBezier

ml = MultiLineString(lines)

pgs = (list(polygonize(ml)))

[plt.plot(*pgs[0].exterior.xy)]
[plt.plot(*ii.xy) for ii in pgs[0].interiors]
plt.show()
