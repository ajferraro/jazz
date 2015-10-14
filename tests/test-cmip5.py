import os

import iris

from jazz import cmip5


def test_path():
    """Test the path generated by cmip5.path exists."""
    path = cmip5.path('HadGEM2-ES', 'historical', 'mon', 'atmos', 'Amon',
                      'r1i1p1', 'ts')
    assert os.path.exists('/'.join(path.split('/')[0:-1]))


def test_fetch():
    path = cmip5.path('IPSL-CM5A-LR', 'amip', 'mon', 'atmos', 'Amon',
                      '*', 'ts')
    cube = cmip5.fetch(path)
    assert isinstance(cube, iris.cube.Cube)


def test_get_fx():
    cube = cmip5.get_fx('IPSL-CM5A-LR', 'sftlf')
    assert isinstance(cube, iris.cube.Cube)
