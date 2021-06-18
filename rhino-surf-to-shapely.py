import rhino3dm as rh
from shapely.geometry import MultiLineString
from shapely.ops import polygonize
import matplotlib.pyplot as plt
import numpy as np
import json
from geomdl import utilities, NURBS, multi
from geomdl.visualization import VisMPL

lines = []

with open('data',"r") as f:
    json_data = f.read()
brep_dict = json.loads(json_data)
brep = rh.CommonObject.Decode(brep_dict[2])

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


mcrv = multi.CurveContainer()
for ii in range(0,len(brep.Edges)):
    nurb = brep.Edges[ii].ToNurbsCurve()
    print(f"ii: {ii}")
    print(f"deg: {nurb.Order-1}")
    print(f"num cp: {len(nurb.Points)}")
    print(f"num knots: {len(nurb.Knots)+2}")
    curve = NURBS.Curve()
    curve.order = nurb.Order
    curve.ctrlpts = [(point.X, point.Y, point.Z) for point in nurb.Points]
    curve.weights = [point.W for point in nurb.Points]
    knots_all = ([nurb.Knots.SuperfluousKnot(True)]+[knot for knot in nurb.Knots]+[nurb.Knots.SuperfluousKnot(False)])
    curve.knotvector = knots_all
    curve.delta = 0.01
    mcrv.add(curve)
vis_config = VisMPL.VisConfig(legend=False, axes=True, figure_dpi=120)
mcrv.vis = VisMPL.VisCurve2D(vis_config)
mcrv.render()
