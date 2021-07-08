from rhino_shapely_interop.importers import RhImporter
import matplotlib.pyplot as plt
from shapely.geometry import Point

rhi = RhImporter.from_file("example_data/file_1.3dm")

# plot all curves on x,y plane
crvs = list(rhi.get_curves(2))
[plt.plot(*ii.xy) for ii in crvs]
plt.show()

# plot all planer brep on x,y plane
polys = list(rhi.get_planer_brep(2))
[plt.plot(*poly.exterior.xy) for poly in polys]
[plt.plot(*interior.xy) for poly in polys for interior in poly.interiors]
plt.show()

# Isolate a single planer brep
pnt = Point(128,40)
sec = None
for poly in polys:
    if poly.contains(pnt): 
        sec = poly
        break

plt.plot(*sec.exterior.xy)
[plt.plot(*interior.xy) for interior in sec.interiors]
plt.show()