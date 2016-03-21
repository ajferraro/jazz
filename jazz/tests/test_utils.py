import numpy as np

import iris
import netcdftime

from .. import utils


def generate_cube():
    """Generates a cube containing month numbers as data, so January always
    has data value 1.

    """
    data = np.array([i+1 for i in xrange(12)]*2)
    dates = np.array([netcdftime.datetime(2000+(i/12), (i % 12)+1, 1)
                      for i in xrange(len(data))])
    time = netcdftime.date2num(dates, 'days since 2000-01-01')
    time_coord = iris.coords.DimCoord(time, standard_name='time', units='days since 2000-01-01')
    dim_coords_and_dims = [(time_coord, 0)]
    cube = iris.cube.Cube(data, dim_coords_and_dims=dim_coords_and_dims,
                          long_name='dummy data')
    return cube


def test_climatology():
    """Climatology should just be month numbers, since the dummy cube is just
    that."""
    cube = generate_cube()
    climatology = utils.climatology(cube)
    ans = np.arange(12)+1
    assert (climatology.data == ans).all()


def test_anomalies():
    """Anomalies should be zero as data are repeated annually."""
    cube = generate_cube()
    anomalies = utils.anomalies(cube)
    assert (anomalies.data == 0).all()
