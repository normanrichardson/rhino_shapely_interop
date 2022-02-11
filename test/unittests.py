import unittest

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

from rhino_shapely_interop.importers import RhImporter
from rhino_shapely_interop.rhino_wrappers import RhCurv, RhPnt
from rhino_shapely_interop.transformations import CoordTransform


class TestCoordTransform(unittest.TestCase):
    def test_simple_cases(self):
        x1 = np.array([1, 0, -1])
        x2 = np.array([0, 1, 0])
        ct = CoordTransform(x1, x2)
        test1 = np.array([[1, 0, 0]])
        test2 = np.array([[0, 1, 0]])
        test3 = np.array([[0, 0, 1]])
        test4 = np.array([[1, 0, 1]])
        test5 = np.array([[0, 1, 1]])
        exp1 = np.array([[0.70710678, 0]])
        exp2 = np.array([[0, 1]])
        exp3 = np.array([[-0.70710678, 0]])
        exp4 = np.array([[0, 0]])
        exp5 = np.array([[-0.70710678, 1]])
        np.testing.assert_array_almost_equal(ct.transform(test1.T), exp1.T)
        np.testing.assert_array_almost_equal(ct.transform(test2.T), exp2.T)
        np.testing.assert_array_almost_equal(ct.transform(test3.T), exp3.T)
        np.testing.assert_array_almost_equal(ct.transform(test4.T), exp4.T)
        np.testing.assert_array_almost_equal(ct.transform(test5.T), exp5.T)

    def test_unit_circle(self):
        n = 101
        theta = np.linspace(0, 2 * np.pi, n)
        pnts = np.array([np.cos(theta), np.sin(theta), np.zeros(n)])
        x1 = np.array([1, 0, -1])
        x2 = np.array([0, 1, 0])
        ct = CoordTransform(x1, x2)
        res = ct.transform(pnts)
        mp = MultiPoint(res.T)
        # test the bounds of the transformed unit circle
        np.testing.assert_array_almost_equal(
            list(mp.bounds), [-1.0 / np.sqrt(2), -1, 1.0 / np.sqrt(2), 1]
        )

    def test_normal(self):
        x1 = np.array([1, 0, 0])
        x2 = np.array([0, 1, 0])
        ct = CoordTransform(x1, x2)
        normal = ct.plane_normal
        exp = np.array([0, 0, 1])
        np.testing.assert_array_almost_equal(normal, exp)


class TestRhCurv(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        uc = rh.ArcCurve(rh.Circle(rh.Point3d(0, 0, 0), 1))
        cls.rc = RhCurv(uc)

    def test_is_line(self):
        self.assertTrue(not self.rc.is_line())

    def test_get_shapely_line(self):
        n = 9
        theta = np.linspace(0, 2 * np.pi, n)
        exp = np.array([np.cos(theta), np.sin(theta)]).T
        test = np.array(self.rc.get_shapely_line(lambda x: x[:2]).coords)
        np.testing.assert_array_almost_equal(test, exp)

    def test_refine(self):
        n = 33
        theta = np.linspace(0, 2 * np.pi, n)
        pnts = np.array([np.cos(theta), np.sin(theta)])
        exp = pnts.T
        self.rc.refine(3)
        test = np.array(
            self.rc.get_shapely_line(transform=lambda x: x[:2]).coords
        )
        np.testing.assert_array_almost_equal(test, exp)

    def tearDown(self):
        self.rc.refine(0)


class TestRhPnt(unittest.TestCase):
    def test_get_shapely_point(self):
        pnt = rh.Point3d(1, 1, 1)
        exp = Point(1, 1)
        rp = RhPnt(pnt)
        test = rp.get_shapely_point(transform=lambda x: x[:2])
        np.testing.assert_array_almost_equal(test.coords, exp.coords)


class TestRhinoImporterValidation(unittest.TestCase):
    def test_validate_file_name(self):
        self.assertTrue(
            RhImporter._validate_file_name("example_data/unit_circle.3dm")
        )

    def test_validate_brep(self):
        pln = rh.Plane(rh.Point3d(0, 0, 0), rh.Vector3d(0, 0, 1))
        bound = rh.Polyline(5)
        bound.Add(1, 1, 0)
        bound.Add(-1, 1, 0)
        bound.Add(-1, -1, 0)
        bound.Add(1, -1, 0)
        bound.Add(1, 1, 0)
        brep = rh.Brep.CreateTrimmedPlane(pln, bound.ToPolylineCurve())
        self.assertTrue(RhImporter._validate_brep(brep))

    def test_validate_curve(self):
        pl = rh.Polyline(3)
        pl.Add(1, 1, 4)
        pl.Add(2, 2, 4)
        pl.Add(3, 1, 4)
        plc = pl.ToPolylineCurve()
        self.assertTrue(RhImporter._validate_curve(plc))

    def test_validate_point(self):
        pnt1 = rh.Point3d(1, 2, 4)
        pnt2 = rh.Point2d(1, 2)
        self.assertTrue(RhImporter._validate_point(pnt1))
        self.assertTrue(RhImporter._validate_point(pnt2))


class TestRhinoImporterCreation(unittest.TestCase):
    def test_creation_from_file(self):
        rhi = RhImporter.from_file("example_data/unit_circle.3dm")
        self.assertEqual(len(rhi._point), 0)
        self.assertEqual(len(rhi._curve), 1)
        self.assertEqual(len(rhi._brep), 0)

    def test_creation_from_file_byte_array(self):
        pass

    def test_creation_from_serialzed_brep(self):
        pln = rh.Plane(rh.Point3d(0, 0, 0), rh.Vector3d(0, 0, 1))
        bound = rh.Polyline(4)
        bound.Add(0, 1, 0)
        bound.Add(-1, -1, 0)
        bound.Add(1, -1, 0)
        bound.Add(0, 1, 0)
        brep = rh.Brep.CreateTrimmedPlane(pln, bound.ToPolylineCurve())
        rhi = RhImporter.from_serialzed_brep(brep.Encode())
        self.assertEqual(len(rhi._point), 0)
        self.assertEqual(len(rhi._curve), 0)
        self.assertEqual(len(rhi._brep), 1)

    def test_creation_from_serialzed_curve(self):
        pl = rh.Polyline(3)
        pl.Add(1, 1, 4)
        pl.Add(2, 2, 4)
        pl.Add(3, 1, 4)
        plc = pl.ToPolylineCurve()
        rhi = RhImporter.from_serialzed_curve(plc.Encode())
        self.assertEqual(len(rhi._point), 0)
        self.assertEqual(len(rhi._curve), 1)
        self.assertEqual(len(rhi._brep), 0)


class TestRhinoImporterExports(unittest.TestCase):
    def test_get_planer_brep(self):
        pln = rh.Plane(rh.Point3d(0, 0, 0), rh.Vector3d(0, 0, 1))
        bound = rh.Polyline(4)
        bound.Add(0, 1, 0)
        bound.Add(-1, -1, 0)
        bound.Add(1, -1, 0)
        bound.Add(0, 1, 0)
        brep = rh.Brep.CreateTrimmedPlane(pln, bound.ToPolylineCurve())
        rhi = RhImporter.from_serialzed_brep(brep.Encode())
        test = [poly for poly in rhi.get_planer_brep(0)][0]
        exp = Polygon([(0, 1), (-1, -1), (1, -1)])
        self.assertTrue((test - exp).is_empty)

    def test_get_curves(self):
        pl = rh.Polyline(5)
        pl.Add(1, 1, 0)
        pl.Add(-1, 1, 0)
        pl.Add(-1, -1, 0)
        pl.Add(1, -1, 0)
        pl.Add(1, 1, 0)
        crv = pl.ToPolylineCurve()
        rhi = RhImporter.from_serialzed_curve(crv.Encode())
        test = [crv for crv in rhi.get_curves(0)][0]
        exp = LineString([(1, 1), (-1, 1), (-1, -1), (1, -1), (1, 1)])
        self.assertTrue((test - exp).is_empty)

    def test_get_points(self):
        pnt1 = rh.Point3d(5, 2, 3)
        pnt2 = rh.Point3d(8, 2, 2)
        model = rh.File3dm()
        model.Objects.AddPoint(pnt1)
        model.Objects.AddPoint(pnt2)
        rhi = RhImporter(model=model)
        test = MultiPoint([poly for poly in rhi.get_points()])
        exp = MultiPoint([(5, 2), (8, 2)])
        self.assertTrue((test - exp).is_empty)
        test = MultiPoint(
            [
                poly
                for poly in rhi.get_points(
                    np.array([0, 0, 1]), np.array([0, 1, 0])
                )
            ]
        )
        exp = MultiPoint([(3, 2), (2, 2)])
        self.assertTrue((test - exp).is_empty)


class TestRhinoImporterLineMergeTol(unittest.TestCase):

    # Test a single line one gap
    def test_single_line_one_gap(self):
        l1 = MultiLineString(
            [
                LineString([(0, 0), (0, 1), (1, 1), (1, 0)]),
                LineString([(1, 0), (0, 0.01)]),
            ]
        )
        l2 = LinearRing([(0, 0.01), (0, 1), (1, 1), (1, 0)])
        l3 = LinearRing([(0, 0), (0, 1), (1, 1), (1, 0)])
        l4 = RhImporter._line_merge_tol(l1, tol=0.02)
        self.assertEqual(len(l4), 1)
        self.assertTrue((l4[0] - l2).is_empty or (l4[0] - l3).is_empty)

    # Test a single line two gaps
    def test_single_line_two_gaps(self):
        l1 = MultiLineString(
            [
                LineString([(0, 0), (0, 1), (1, 1), (1, 0)]),
                LineString([(1.01, 0), (0, 0.01)]),
            ]
        )
        l2 = LinearRing([(0, 0.01), (0, 1), (1, 1), (1.01, 0)])
        l3 = LinearRing([(0, 0), (0, 1), (1, 1), (1, 0)])
        l4 = RhImporter._line_merge_tol(l1, tol=0.02)
        self.assertEqual(len(l4), 1)
        self.assertTrue((l4[0] - l2).is_empty or (l4[0] - l3).is_empty)

    # Test multiple lines one gap each
    def test_multiple_line_one_gap(self):
        l1 = MultiLineString(
            [
                LineString([(-1, -1), (-1, 1), (1, 1), (1, -1)]),
                LineString([(1, -1), (-1, -1.01)]),
                LineString([(-2, -2), (-2, 2), (2, 2), (2, -2)]),
                LineString([(2, -2), (-2, -2.01)]),
            ]
        )
        l2 = LinearRing([(-1, -1), (-1, 1), (1, 1), (1, -1)])
        l3 = LinearRing([(-2, -2), (-2, 2), (2, 2), (2, -2)])
        l4 = RhImporter._line_merge_tol(l1, tol=0.02)
        self.assertEqual(len(l4), 2)
        np.testing.assert_array_almost_equal(
            np.asarray(l2.coords), np.asarray(l4[1].coords), 2
        )
        np.testing.assert_array_almost_equal(
            np.asarray(l3.coords), np.asarray(l4[0].coords), 2
        )

    # Test multiple lines two gap each
    def test_multiple_line_two_gap(self):
        gap = 0.01
        l1 = MultiLineString(
            [
                LineString([(-1, -1), (-1, 1), (1, 1), (1, -1-gap)]),
                LineString([(1, -1), (-1, -1+gap)]),
                LineString([(-2, -2), (-2, 2), (2, 2), (2, -2-gap)]),
                LineString([(2, -2), (-2, -2+gap)]),
            ]
        )
        eps = np.finfo(float).eps
        l2 = LinearRing([(-1, -1), (-1, 1), (1, 1), (1, -1)])
        l3 = LinearRing([(-2, -2), (-2, 2), (2, 2), (2, -2)])
        l4 = RhImporter._line_merge_tol(l1, tol=gap+eps)
        self.assertEqual(len(l4), 2)
        np.testing.assert_array_almost_equal(
            np.asarray(l2.coords), np.asarray(l4[1].coords), 2
        )
        np.testing.assert_array_almost_equal(
            np.asarray(l3.coords), np.asarray(l4[0].coords), 2
        )

    # Test tol too large causing endpoints to desolve
    def test_tol_too_large(self):
        gap = 0.01
        eps = np.finfo(float).eps
        p1 = Point(-1, -1+gap)
        p2 = Point(-2, -2+gap)
        tol = p1.distance(p2)
        tol_big = p1.distance(p2)+eps

        l1 = MultiLineString(
            [
                LineString([(-1, -1), (-1, 1), (1, 1), (1, -1)]),
                LineString([(1, -1), p1]),
                LineString([(-2, -2), (-2, 2), (2, 2), (2, -2)]),
                LineString([(2, -2), p2]),
            ]
        )

        # Test that the expected tolerance works
        test = False
        try:
            RhImporter._line_merge_tol(l1, tol=tol)
        except ValueError:
            test = True
        self.assertFalse(test, 'Expected accepted tolerance failed')

        # Test that the expected failed tolerance (eps larger than then previous tolerance) fails 
        with self.assertRaises(ValueError) as cm:
            RhImporter._line_merge_tol(l1, tol=tol_big)
        self.assertEqual(str(cm.exception), "Multiple endpints snaped. Individual lines are 'desolving'. Tolerance is too high.")

        # Test tol too large causing endpoints to desolve
    def test_tol_too_small(self):
        gap = 0.01
        eps = np.finfo(float).eps
        p1 = Point(-1, -1+gap)
        p2 = Point(-2, -2+gap)
        tol = gap+eps
        tol_small = gap

        l1 = MultiLineString(
            [
                LineString([(-1, -1), (-1, 1), (1, 1), (1, -1)]),
                LineString([(1, -1), p1]),
                LineString([(-2, -2), (-2, 2), (2, 2), (2, -2)]),
                LineString([(2, -2), p2]),
            ]
        )

        # Test that the expected tolerance works
        test = False
        try:
            RhImporter._line_merge_tol(l1, tol=tol)
        except ValueError:
            test = True
        self.assertFalse(test, 'Expected accepted tolerance failed')

        # Test that the expected failed tolerance (eps larger than then previous tolerance) fails 
        with self.assertRaises(ValueError) as cm:
            RhImporter._line_merge_tol(l1, tol=tol_small)
        self.assertEqual(str(cm.exception), "The manual line merged failed to form LinearRings within the allowable tolerance. Tolerance too small.")

if __name__ == "__main__":
    unittest.main()
