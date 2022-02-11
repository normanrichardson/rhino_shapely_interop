import os.path
import uuid
from typing import Iterator, List, Optional, Union

import numpy as np
import rhino3dm as rh
from shapely.geometry import (
    LinearRing,
    LineString,
    MultiLineString,
    MultiPoint,
    Point,
    Polygon,
)
from shapely.ops import linemerge, polygonize, snap

from .rhino_wrappers import RhCurv, RhPnt
from .transformations import CoordTransform


class RhImporter:
    """Import geometric objects from rhino.

    Parameters
    ----------
    model : rhino3dm.File3dm
        A keyword argument for a rhino3dm model
    brep : rhino3dm.Brep
        A keyword argument for a rhino3dm brep
    curve : rhino3dm.Curve
        A keyword argument for a rhino3dm curve


    Methods
    -------
    from_file(file_name) :
        Classmethod to create the object from a file.
    from_file_byte_array() :
        Classmethod to create the object from a byte array file object.
    from_serialzed_brep(s_brep) :
        Classmethod to create the object from a serialized brep object.
    from_serialzed_curve(s_curve) :
        Classmethod to create the object from a serialized curve object.
    get_planer_brep(refine_num, vec1, vec2, plane_distance, project, parallel) :
        Generator that returns single surface planer breps as shapely polygons.
    get_curves(refine_num, vec1, vec2, plane_distance, project, parallel) :
        Generator that returns curves as shapely line strings.
    get_points(vec1, vec2, plane_distance, project) :
        Generator that returns points as shapely points.
    """

    def __init__(
        self,
        *,
        model: Optional[rh.File3dm] = None,
        brep: Optional[rh.Brep] = None,
        curve: Optional[rh.Curve] = None
    ):
        """Constructor

        Parameters
        ----------
        model : rhino3dm.File3dm
            A keyword argument for a rhino3dm model
        brep : rhino3dm.Brep
            A keyword argument for a rhino3dm brep
        curve : rhino3dm.Curve
            A keyword argument for a rhino3dm curve
        """
        self._brep = []
        self._curve = []
        self._point = []
        if model is not None:
            inst_ids = self._get_instance_ids(model)
            self._process_objects(model.Objects, inst_ids)
        if brep is not None:
            self._brep.append(brep)
        if curve is not None:
            self._curve.append(curve)

    @classmethod
    def from_file(cls, file_name: str) -> "RhImporter":
        """Class method to read from a Rhino file.

        Parameters
        ----------
        file_name : string
            File name of (or path to) a Rhino file.

        Returns
        -------
        RhImporter

        Raises
        ------
        ValueError
            If file extention is not ".3dm" or the file does not exit.
        """
        if not cls._validate_file_name(file_name):
            raise ValueError("File does not exist or is not a Rhino file.")
        model = rh.File3dm.Read(file_name)
        return cls(model=model)

    @classmethod
    def from_file_byte_array(cls, s_file: str) -> "RhImporter":
        """Class method to read from a Rhino file from a byte array.

        Parameters
        ----------
        file_name : bytes
            Byte array of file

        Returns
        -------
        RhImporter

        Raises
        ------
        ValueError
            If the byte array could not be serialized.
        """
        model = rh.File3dm.FromByteArray(s_file)
        if model is None:
            ValueError("Byte array could not be serialized.")
        return cls(model=model)

    @classmethod
    def from_serialzed_brep(cls, s_brep: str) -> "RhImporter":
        """Class method to read from a serialized brep object.

        Parameters
        ----------
        s_brep : string
            Serialization of a brep

        Returns
        -------
        RhImporter

        Raises
        ------
        ValueError
            If the provided serialized object is not a single surface planer
            brep
        """
        brep = rh.CommonObject.Decode(s_brep)
        if not cls._validate_brep(brep):
            raise ValueError("Data is not a single surface planer brep.")
        return cls(brep=brep)

    @classmethod
    def from_serialzed_curve(cls, s_curve: str) -> "RhImporter":
        """Class method to read from a serialized curve object.

        Parameters
        ----------
        s_curve : string
            Serialization of a curve

        Returns
        -------
        RhImporter

        Raises
        ------
        ValueError
            If the provided serialized object is not a curve
        """
        curve = rh.CommonObject.Decode(s_curve)
        if not cls._validate_curve(curve):
            raise ValueError("Data is not a curve.")
        return cls(curve=curve)

    @staticmethod
    def _validate_file_name(file_name: str) -> bool:
        """Performs checks on a file name.

        Parameters
        ----------
        file_name : string
            The file name

        Returns
        -------
        Boolean
        """
        valid = False
        if os.path.isfile(file_name):
            valid = file_name.endswith(".3dm")
        return valid

    @staticmethod
    def _validate_brep(geom: rh.Brep) -> bool:
        """Performs checks on a brep.

        Parameters
        ----------
        geom : rhino3dm.GeometryBase
            Geometry to test for brep properties.

        Returns
        -------
        Boolean
        """
        valid = False
        if isinstance(geom, rh.Brep):
            if len(geom.Surfaces) == 1:
                valid = geom.Surfaces[0].IsPlanar()
        return valid

    @staticmethod
    def _validate_curve(geom: rh.Curve) -> bool:
        """Performs checks on a curve.

        Parameters
        ----------
        geom : rhino3dm.GeometryBase
            Geometry to test for curve properties.

        Returns
        -------
        Boolean
        """
        return isinstance(geom, rh.Curve)

    @staticmethod
    def _validate_point(geom: rh.Point) -> bool:
        """Performs checks on a point.

        Parameters
        ----------
        geom : rhino3dm.GeometryBase
            Geometry to test for point properties.

        Returns
        -------
        Boolean
        """
        return isinstance(geom, (rh.Point3d, rh.Point, rh.Point2d))

    @staticmethod
    def _get_instance_ids(model: rh.File3dm) -> List[uuid.UUID]:
        """
        Gets the list of geometric objects that are part of the models
        InstanceDefinitions. These are geometric objects that are:
        * Blocks
        * Annotations

        Parameters
        ----------
        model : rhino3dm.File3dm
            A rhino3dm model

        Returns
        -------
        List of UUID
            List of geometric objects that are part of the models
            InstanceDefinitions
        """
        ids = []
        for instance in model.InstanceDefinitions:
            ids += instance.GetObjectIds()
        return ids

    def _process_objects(
        self, objects: rh.File3dmObjectTable, inst_ids: List[uuid.UUID]
    ) -> None:
        """Process the file object table.

        Parameters
        ----------
        objects : rhino3dm.File3dmObjectTable
            The object table for the file.
        """
        for obj in objects:
            if obj.Attributes.Id not in inst_ids:
                if self._validate_brep(obj.Geometry):
                    self._brep.append(obj.Geometry)
                elif self._validate_curve(obj.Geometry):
                    self._curve.append(obj.Geometry)
                elif self._validate_point(obj.Geometry):
                    self._point.append(obj.Geometry)

    def get_planer_brep(
        self,
        refine_num: Optional[int] = 1,
        tol: float = 1e-10,
        vec1: Optional[np.ndarray] = np.array([1, 0, 0]),
        vec2: Optional[np.ndarray] = np.array([0, 1, 0]),
        plane_distance: Optional[float] = 0.0,
        project: Optional[bool] = True,
        parallel: Optional[bool] = False,
    ) -> Iterator[Polygon]:
        """Get all the single surface planer breps as Shapely polygons.
        Two vectors `vec1` and `vec2` describe the Shapely plane, with
        coordinates (x',y'). The breps coordinates (x,y,z) are projected onto
        (x',y'). Options to filter breps are provided:
        * only breps that are parallel to the Shapely plane
        * only breps in the Shapely plane

        Parameters
        ----------
        refine_num : integer, optional
            Bézier curve interpolation number. In Rhino a surface's edges are
            nurb based curves. Shapely does not support nurbs, so the
            individual Bézier curves are interpolated using straight lines.
            This parameter sets the number of straight lines used in the
            interpolation.
        tol : float, optional
            A tolerance used to manually merge the endpoints of brep edges.
            This only occurs if the brep edges do not merge to LinearRings.
            Consider using the model's absolute tolerance:
            rhino3dm.File3dmSettings.ModelAbsoluteTolerance.
        vec1 : numpy array, optional
            A 3d vector in the Shapely plane. Rhino is a 3D geometry
            environment. Shapely is a 2D geometric library. Thus a 2D plane
            needs to be defined in Rhino that represents the Shapely coordinate
            system. `vec1` represents the 1st vector of this plane. It will be
            used as Shapely's x direction.
        vec2 : numpy array, optional
            Continuing from `vec1`, `vec2` is another vector to define the
            Shapely plane. It must not be [0,0,0] and it's only requirement is
            that it is any vector in the Shapely plane (but not equal to
            `vec1`).
        plane_distance : float, optional
            The distance to the Shapely plane.
        project : Boolean, optional
            Controls if the breps are projected onto the plane in the direction
            of the Shapley plane's normal.
        parallel : Boolean, optional
            Controls if only the rhino surfaces that have the same normal as
            the Shapely plane are yielded. If true, all non parallel surfaces
            are filtered out.

        Yields
        -------
        Shapely polygon.
            Shapely polygons with a coordinate system defined by vec1, and
            vec2.

        Raises
        ------
        ValueError
            If vec1 is not of ndim=1 and size=3
        ValueError
            If vec2 is not of ndim=1 and size=3
        ValueError
            If vec1 == vec2
        ValueError
            If vec1 or vec2 are the origin
        ValueError
            If both project and parallel are false.
        """

        if not (vec1.ndim == 1 and vec1.size == 3):
            raise ValueError("vec1 is a numpy vector in 3d")
        if not (vec2.ndim == 1 and vec2.size == 3):
            raise ValueError("vec2 is a numpy vector in 3d")
        if not project and not parallel:
            raise ValueError(
                (
                    "No surface meets this criteria, a surface that is not "
                    "parallel and is not projected. This would just be the "
                    "intersction of the plane and the surface (i.e. a line)."
                )
            )
        if (vec1 == vec2).all():
            raise ValueError("vec2 must be different from vec1.")
        if (vec1 == 0).all() or (vec2 == 0).all():
            raise ValueError("The vectors must not be the origin.")

        ct = CoordTransform(vec1, vec2)

        def validation_factory():
            if project and parallel:
                return lambda normal, *args: (
                    ct.plane_normal == np.array([normal.X, normal.Y, normal.Z])
                ).all()
            if project and not parallel:
                return lambda *args: True
            if not project and parallel:
                return (
                    lambda normal, orig: (
                        ct.plane_normal
                        == np.array([normal.X, normal.Y, normal.Z])
                    ).all()
                    and normal.X * orig.X
                    + normal.Y * orig.Y
                    + normal.Z * orig.Z
                    == plane_distance
                )

        validation = validation_factory()

        for brep in self._brep:
            _, frame = brep.Surfaces[0].FrameAt(0, 0)
            if validation(frame.ZAxis, frame.Origin):
                rh_curvs = []
                for ii in range(0, len(brep.Edges)):
                    rh_curvs.append(RhCurv(brep.Edges[ii]))
                for rh_curv in rh_curvs:
                    if not rh_curv.is_line():
                        rh_curv.refine(refine_num)

                pw_line_list = MultiLineString(
                    [rc.get_shapely_line(ct.transform) for rc in rh_curvs]
                )
                line_list = linemerge(pw_line_list)

                # When deconstructing the brep into individual edges the
                # endpoints of the curves are subject to a precision. This
                # leads to Shapely's linemerge algorithm to fail to create
                # LinearRings, as there are small gaps in the lines.
                line_list = self._line_merge_tol(line_list, tol)

                pgs = list(polygonize(line_list))

                if len(pgs) == 1:
                    yield pgs[0]
                # filter out the holes returned in polygonize. There should be
                # only one ploygon returned that has interiors
                if len(pgs) > 0:
                    for poly in pgs:
                        if len(poly.interiors) > 0:
                            yield poly

    def get_curves(
        self,
        refine_num: Optional[int] = 1,
        vec1: Optional[np.ndarray] = np.array([1, 0, 0]),
        vec2: Optional[np.ndarray] = np.array([0, 1, 0]),
        plane_distance: Optional[float] = 0.0,
        project: Optional[bool] = True,
        parallel: Optional[bool] = False,
    ) -> Iterator[LineString]:
        """Get all rhino curves as Shapely line strings.
        Two vectors `vec1` and `vec2` describe the Shapely plane, with
        coordinates (x',y'). The rhino curves coordinates (x,y,z) are
        projected onto (x',y'). Options to filter curve are provided:
        * only curves that are parallel to the Shapely plane
        * only curves in the Shapely plane

        Parameters
        ----------
        refine_num : integer
            Bézier curve interpolation number. In Rhino a curves are nurb
            based. Shapely does not support nurbs, so the individual Bézier
            curves are interpolated using straight lines. This parameter sets
            the number of straight lines used in the interpolation.
        vec1 : numpy array, optional
            A 3d vector in the Shapely plane. Rhino is a 3D geometry
            environment. Shapely is a 2D geometric library. Thus a 2D plane
            needs to be defined in Rhino that represents the Shapely coordinate
            system.
            `vec1` represents the 1st vector of this plane. It will be used as
            Shapely's x direction.
        vec2 : numpy array, optional
            Continuing from `vec1`, `vec2` is another vector to define the
            Shapely plane. It must not be [0,0,0] and it's only requirement is
            that it is any vector in the Shapely plane (but not equal to
            `vec1`).
        plane_distance : float, optional
            The distance to the Shapely plane.
        project : Boolean, optional
            Controls if the rhino curves are projected onto the plane in the
            direction of the Shapley plane's normal.
        parallel : Boolean, optional
            Controls if the yield curves are in a plane parallel to the Shapely
            plane. If true, all non parallel surfaces are filtered out.

        Yields
        -------
        Shapely lineString.
            Shapely lineString with a coordinate system defined by vec1, and
            vec2.

        Raises
        ------
        ValueError
            If vec1 is not of ndim=1 and size=3
        ValueError
            If vec2 is not of ndim=1 and size=3
        ValueError
            If vec1 == vec2
        ValueError
            If vec1 or vec2 are the origin
        ValueError
            If both project and parallel are false.
        """

        if not (vec1.ndim == 1 and vec1.size == 3):
            raise ValueError("vec1 is a numpy vector in 3d")
        if not (vec2.ndim == 1 and vec2.size == 3):
            raise ValueError("vec2 is a numpy vector in 3d")
        if not project and not parallel:
            raise ValueError(
                (
                    "No surface meets this criteria, a surface that is not "
                    "parallel and is not projected. This would just be the "
                    "intersction of the plane and the surface (i.e. a line)."
                )
            )
        if (vec1 == vec2).all():
            raise ValueError("vec2 must be different from vec1.")
        if (vec1 == 0).all() or (vec2 == 0).all():
            raise ValueError("The vectors must not be the origin.")

        ct = CoordTransform(vec1, vec2)

        def validation_factory():
            if project and parallel:
                # 1)check it is planer
                # 2) check 2 points have the same distance value
                return (
                    lambda planer, *args: planer
                    and ct.plane_normal.dot(args[0])
                    - ct.plane_normal.dot(args[1])
                    == 0
                )
            if project and not parallel:
                return lambda *args: True
            if not project and parallel:
                # 1)check it is planer 2) check a point is in the plane
                return (
                    lambda planer, *args: planer
                    and ct.plane_normal.dot(args[0]) == plane_distance
                )

        validation = validation_factory()

        for curve in self._curve:
            curve_w = RhCurv(curve)
            if validation(curve_w.is_planer, *curve_w.get_greville_points):
                if not curve_w.is_line():
                    curve_w.refine(refine_num)
                ls = curve_w.get_shapely_line(ct.transform)
                yield ls

    def get_points(
        self,
        vec1: Optional[np.ndarray] = np.array([1, 0, 0]),
        vec2: Optional[np.ndarray] = np.array([0, 1, 0]),
        plane_distance: Optional[float] = 0.0,
        project: Optional[bool] = True,
    ) -> Iterator[Point]:
        """Get all the rhino points as Shapely points.
        Two vectors `vec1` and `vec2` describe the Shapely plane, with
        coordinates (x',y'). The point coordinates (x,y,z) are projected onto
        (x',y'). Options to filter rhino points are provided:
        * only points in the Shapely plane

        Parameters
        ----------
        vec1 : numpy array, optional
            A 3d vector in the Shapely plane. Rhino is a 3D geometry
            environment. Shapely is a 2D geometric library.
            Thus a 2D plane needs to be defined in Rhino that represents the
            Shapely coordinate system. `vec1` represents the 1st vector of this
            plane. It will be used as Shapely's x direction.
        vec2 : numpy array, optional
            Continuing from `vec1`, `vec2` is another vector to define the
            Shapely plane. It must not be [0,0,0] and it's only requirement is
            that it is any vector in the Shapely plane (but not equal to
            `vec1`).
        plane_distance : float, optional
            The distance to the Shapely plane.
        project : Boolean, optional
            Controls if the points are projected onto the plane in the
            direction of the Shapley plane's normal.

        Yields
        -------
        Shapely point.
            Shapely points with a coordinate system defined by vec1, and vec2.

        Raises
        ------
        ValueError
            If vec1 is not of ndim=1 and size=3
        ValueError
            If vec2 is not of ndim=1 and size=3
        ValueError
            If vec1 == vec2
        ValueError
            If vec1 or vec2 are the origin
        """

        if not (vec1.ndim == 1 and vec1.size == 3):
            raise ValueError("vec1 is a numpy vector in 3d")
        if not (vec2.ndim == 1 and vec2.size == 3):
            raise ValueError("vec2 is a numpy vector in 3d")
        if (vec1 == vec2).all():
            raise ValueError("vec2 must be different from vec1.")
        if (vec1 == 0).all() or (vec2 == 0).all():
            raise ValueError("The vectors must not be the origin.")

        ct = CoordTransform(vec1, vec2)

        def validation_factory():
            if project:
                return lambda *args: True
            if not project:
                return lambda pnt: ct.plane_normal.dot(pnt) == -plane_distance

        validation = validation_factory()
        for pnt in self._point:
            pnt_w = RhPnt(pnt)
            if validation(pnt_w.as_numpy):
                yield pnt_w.get_shapely_point(ct.transform)

    @staticmethod
    def _line_merge_tol(
        crvs: Union[LineString, MultiLineString], tol: float = 1e-12
    ) -> MultiLineString:

        completed = []
        incomplete = []
        # Search for rings
        try:
            # If a MultiLineString
            for crv in crvs.geoms:
                if crv.is_ring:
                    completed.append(crv)
                else:
                    incomplete.append(crv)
        except AttributeError:
            # If a LineString
            if crvs.is_ring:
                completed.append(crvs)
            else:
                completed.append(LinearRing(crvs))

        # If all curvs are rings or only one LineSting, exit
        if len(incomplete) == 0:
            return completed

        # Try snap the end points of lines in the incomplete group to one
        # another end points of lines (snap the end points based on a
        # tollerance)
        processed = []
        for (idx, crv_i) in enumerate(incomplete):
            comp = crv_i
            snaps_a = 0
            snaps_b = 0
            for crv_j in incomplete[idx + 1 :]:
                p_a = Point(comp.coords[0])
                p_b = Point(comp.coords[-1])
                p_targets = MultiPoint([crv_j.coords[0], crv_j.coords[-1]])
                p_a_new = snap(p_a, p_targets, tolerance=tol)
                p_b_new = snap(p_b, p_targets, tolerance=tol)
                # if a snap happened count it
                if not (p_a - p_a_new).is_empty:
                    snaps_a += 1
                if not (p_b - p_b_new).is_empty:
                    snaps_b += 1
                # if more than one snap occurs raise an error
                if snaps_a > 1 or snaps_b > 1:
                    raise ValueError(
                        "Multiple endpints snaped. Individual lines are "
                        "'desolving'. Tolerance is too high."
                    )
                p_a, p_b = p_a_new, p_b_new
                comp = LineString([p_a, *comp.coords[1:-1], p_b])
            processed.append(comp)

        # Merge the snapped lines
        merge = linemerge(processed)

        # Search for rings in the merged lines
        try:
            # If a MultiLineString
            for crv in merge.geoms:
                if crv.is_ring:
                    completed.append(crv)
                else:
                    d = Point(crv.coords[0]).distance(Point(crv.coords[-1]))
                    if d > tol:
                        raise ValueError(
                            "The manual line merged failed to form LinearRings "
                            "within the allowable tolerance. "
                            "Tolerance too small."
                        )
                    completed.append(LinearRing(crv))
        except AttributeError:
            # If a LineString
            if merge.is_ring:
                completed.append(merge)
            else:
                d = Point(merge.coords[0]).distance(Point(merge.coords[-1]))
                if d > tol:
                    raise ValueError(
                        "The manual line merged failed to form LinearRings "
                        "within the allowable tolerance. Tolerance too small."
                    )
                completed.append(LinearRing(merge))

        return completed
