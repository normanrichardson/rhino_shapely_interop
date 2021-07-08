# rhino_shapely_interop
A python package for converting rhino geometry (https://www.rhino3d.com/) into shapely geometry (https://pypi.org/project/Shapely/).

## Supported Geometry:
1) Rhino Points -> Shapely Points
2) Rhino Curves -> Shapely LineStrings
3) Rhino Breps  -> Shapely Polygon (limited to planer single surface breps)

## Current Limitations:
1) One way (Rhino to Shapely)
2) Rhino is 3D, Shapely is 2D (effectively)
    * Breps are limited to planer single surface breps
    * A shapely plane is defined 
3) InstanceDefinition geometries are ignored (geometry listed in rhino3dm.File3dmInstanceDefinitionTable)
    * Blocks
    * Annotations

## Testing:
python -m unittest test.unittests

## Examples:
See 
1) `example1.py` for importing from a file
2) `example2.py` for importing from serialized brep