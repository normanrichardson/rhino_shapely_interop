Change Log:
==========

v0.0.4:
-------

- Extend the algorithm in RhImporter.get_planer_brep. 
  Issue `#7 <https://github.com/normanrichardson/rhino_shapely_interop/pull/7>`_ illustrated that the shapely linemerge can fail due to precision error. These errors are most likely caused by floating-point representation on rhino3dm.

  For this fix `#8 <https://github.com/normanrichardson/rhino_shapely_interop/pull/8>`_, additional steps are introduced after linemerge. These steps try to snap the endpoints of the edges together. The snapping requires a provided tolerance. The routine can raise errors if the tolerance is too large or too small. A large tolerance will mean that the progressive snapping 'dissolves' the endpoints. A small tolerance will result in the brep edges not resolving to LinearRings or LineStrings (whose endpoints are within the tolerance).
  
  A new optional argument has been introduced to RhImporter.get_planer_brep, tol for the snapping tolerance.

v0.0.3:
-------

- Updates assocciated with changes in Shapely 1.8.0.

v0.0.2:
-------

- Formalize dependency requirements. Mainly slacken constraints on distribution dependencies.

v0.0.1:
-------

- Initial release of the package.
