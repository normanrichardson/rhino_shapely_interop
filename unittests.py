from shapely.geometry.multipoint import asMultiPoint
from rhino_to_shapely import CoordTransform
import numpy as np
import unittest
from shapely.geometry import Point, asMultiPoint, MultiPoint
from shapely.affinity import affine_transform
import matplotlib.pyplot as plt

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
        
if __name__=='__main__':
      unittest.main()