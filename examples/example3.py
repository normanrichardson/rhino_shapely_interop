import matplotlib.pyplot as plt

from rhino_shapely_interop.importers import RhImporter

rhi = RhImporter.from_file("example_data/DIES01_POLYLINES.3dm")

# plot all planer brep on x,y plane
polys = list(rhi.get_planer_brep(2, tol=10))
[plt.plot(*poly.exterior.xy) for poly in polys]
[plt.plot(*interior.xy) for poly in polys for interior in poly.interiors]
plt.show()
