# rhino_shapely_interop
A python package for converting rhino geometry (https://www.rhino3d.com/) into shapely geometry (https://pypi.org/project/Shapely/).

## Installation
> pip install rhino-shapley-interop

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

## Examples:
See 
1) [`example1.py`](https://github.com/normanrichardson/rhino_shapely_interop/blob/master/example1.py) for importing from a file
2) [`example2.py`](https://github.com/normanrichardson/rhino_shapely_interop/blob/master/example2.py) for importing from serialized brep

# Development/Contributions
1. Fork and clone to a local working directory
2. Setup a virtual environment
> python -m venv env

> source env/bin/activate
3. Install in editable mode
> pip install -e .

4. Testing
> python -m unittest test.unittests