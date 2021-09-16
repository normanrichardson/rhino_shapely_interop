from rhino_shapely_interop.importers import RhImporter
import matplotlib.pyplot as plt
import json

# load the encoded json object of a brep object
with open("example_data/data.json") as file:
    data = json.load(file)

rhi = RhImporter.from_serialzed_brep(data)

# plot all planer brep on x,y plane
polys = list(rhi.get_planer_brep(2))
[plt.plot(*poly.exterior.xy) for poly in polys]
[plt.plot(*interior.xy) for poly in polys for interior in poly.interiors]
plt.show()