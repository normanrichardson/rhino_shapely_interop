import numpy as np


class CoordTransform:
    """Class that is responsible for the tranformations between a
    3d environment (x,y,z) and a 2d environment (x',y'). The 2d plane is
    derived from two vectors that define the 2d plane.

    Parameters
    ----------
    vec1 : ndarray
        A 3d vector (x,y,z) in the prefered 2d plane. This is normalised and
        used as e1.
    vec2 : ndarray
        Another 3d vector (x,y,z) in the prefered 2d plane.

    Attributes
    ----------
    plane_normal : ndarray
        The normal of the plane in (x,y,z).

    Methods
    -------
    transform(pnts) :
        Transforms a coordinate from the (x,y,z) into (x',y')
    plane_normal() :
        The vector normal to `vec1` and `vec2`.
    """
    def __init__(self, vec1: np.ndarray, vec2: np.ndarray):
        """Constructor

        Parameters
        ----------
        vec1 : ndarray
            A 3d vector (x,y,z) in the prefered 2d plane. This is normalised
            and used as e1.
        vec2 : ndarray
            Another 3d vector (x,y,z) in the prefered 2d plane.
        """
        self._e1 = vec1 / np.linalg.norm(vec1)
        e3 = np.cross(vec1, vec2)
        self._e3 = e3 / np.linalg.norm(e3)
        self._e2 = np.cross(self._e3, self._e1)
        self.Tinv = np.array([self._e1, self._e2, self._e3]).T
        self.T = np.linalg.inv(self.Tinv)

    def transform(self, pnts: np.ndarray) -> np.ndarray:
        """Transforms a coordinate from the (x,y,z) into (x',y')

        Parameters
        ----------
        pnts : ndarray
            Points in (x,y,z)

        Returns
        -------
        ndarray
            Points in (x',y')
        """
        pnts_prime = self.T.dot(pnts)
        if pnts.ndim == 1:
            return pnts_prime[0:2]
        return pnts_prime[0:2, :]

    @property
    def plane_normal(self) -> np.ndarray:
        """The vector normal to `vec1` and `vec2`.

        Returns
        -------
        ndarray
        """
        return self._e3
