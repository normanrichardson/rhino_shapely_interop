from rhino_to_shapely import CoordTransform, RhCurv
import unittest
from shapely.geometry import Point, asMultiPoint
import numpy as np
import rhino3dm as rh

class TestCoordTransform(unittest.TestCase):
    
    def test_simple_cases(self):
        x1 = np.array([1,0,-1])
        x2 = np.array([0,1,0])
        ct = CoordTransform(x1,x2)
        test1 = np.array([[1,0,0]])
        test2 = np.array([[0,1,0]])
        test3 = np.array([[0,0,1]])
        test4 = np.array([[1,0,1]])
        test5 = np.array([[0,1,1]])
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
        theta = np.linspace(0,2*np.pi, n)
        pnts = np.array([np.cos(theta), np.sin(theta), np.zeros(n)])
        x1 = np.array([1,0,-1])
        x2 = np.array([0,1,0])
        ct = CoordTransform(x1,x2)
        res = ct.transform(pnts)
        mp = asMultiPoint(res.T)
        # test the bounds of the transformed unit circle
        np.testing.assert_array_almost_equal(list(mp.bounds), [-1.0/np.sqrt(2), -1, 1.0/np.sqrt(2), 1])
    
    def test_normal(self):
        x1 = np.array([1,0,0])
        x2 = np.array([0,1,0])
        ct = CoordTransform(x1,x2)
        normal = ct.plane_normal
        exp = np.array([0,0, 1])
        np.testing.assert_array_almost_equal(normal, exp)

class TestRhCurv(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        uc = rh.ArcCurve(rh.Circle(rh.Point3d(0,0,0), 1))
        cls.rc = RhCurv(uc)

    def test_is_line(self):
        self.assertTrue(not self.rc.is_line())

    def test_get_shapely_line(self):
        n = 9
        theta = np.linspace(0,2*np.pi, n)
        exp = np.array([np.cos(theta), np.sin(theta)]).T
        test = np.array(self.rc.get_shapely_line(lambda x: x[:2]))
        np.testing.assert_array_almost_equal(test, exp)

    def test_refine(self):
        n = 33
        theta = np.linspace(0,2*np.pi, n)
        pnts = np.array([np.cos(theta), np.sin(theta)])
        exp = pnts.T
        self.rc.refine(3)
        test = np.array(self.rc.get_shapely_line(lambda x: x[:2]))
        np.testing.assert_array_almost_equal(test, exp)

    def tearDown(self):
        self.rc.refine(0)

if __name__=='__main__':
      unittest.main()