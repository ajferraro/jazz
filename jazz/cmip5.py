"""
Module for handling CMIP5 files
"""
import glob
import itertools
import numpy as np
import os
import warnings

import iris
import iris.coords

import names
import utils


def path(model, experiment, time_frequency, realm, cmor_table, ensemble,
         variable):
    """Return wild-card path for CMIP5 files.

    Args:
        experiment (str): CMOR experiment ID
        model (str): CMOR model ID
        variable (str): CMOR variable ID
        ensemble (str, optional): ensemble ID

    Returns:
        str

    """
    geomip_experiments = ['G1', 'G1oceanAlbedo', 'G2', 'G3', 'G3S', 'G4',
                          'G4cdnc', 'G4seaSalt']
    if experiment in geomip_experiments:
        cmip5_dir = '/badc/cmip5/data/GeoMIP/output1/'
    else:
        cmip5_dir = '/badc/cmip5/data/cmip5/output1/'

    from names import institutes
    institute = institutes[model]

    path = os.path.join(cmip5_dir, institute, model, experiment,
                        time_frequency, realm, cmor_table, ensemble, 'latest',
                        variable, '*')
    return path


def get_fx(model, variable, constraint=None):
    """Get data from the fx table.  Data from this table are usually the same
    in all simulatins (orography, land fraction, etc.).  They are therefore
    not always stored in the directories for each experiment.  This function
    looks in several places to find the data.

    Args:
        model (str): CMIP5 model name
        variable (str): CMIP5 variable name (from the fx CMOR table)

    Returns:
        iris.cube.Cube

    """
    try:
        data = fetch(path(names.coupled_model(model), 'historical', 'fx',
                          'atmos', 'fx', 'r0i0p0', variable), constraint)
    except:
        try:
            data = fetch(path(names.coupled_model(model), 'piControl', 'fx',
                              'atmos', 'fx', 'r0i0p0', variable), constraint)
        except:
            data = fetch(path(names.atmos_model(model), 'amip', 'fx',
                              'atmos', 'fx', 'r0i0p0', variable), constraint)
    return data


def get_orog(model):
    """Return path to access orography data for CMIP5 model.

    This needs a special function because some the files can be found in
    either historical or piControl directories.

    Args:
        model (str): CMIP5 model name

    Returns:
        str

    """

    if model == 'CanAM4':
        coupled_model = 'CanESM2'
    elif model == 'HadGEM2-A':
        coupled_model = 'HadGEM2-ES'
    else:
        coupled_model = model

    try:
        orog = fetch(path(coupled_model, 'historical', 'fx', 'atmos', 'fx',
                          'r0i0p0', 'orog'))
    except:
        try:
            orog = fetch(path(coupled_model, 'piControl', 'fx', 'atmos', 'fx',
                              'r0i0p0', 'orog'))
        except:
            orog = fetch(path(model, 'amip', 'fx', 'atmos', 'fx',
                              'r0i0p0', 'orog'))
    return orog


def availability(model, experiment, time_frequency, realm, cmor_table,
                 ensemble, variable):
    return sorted(glob.glob('/badc/cmip5/data/cmip5/output1/*/{0}/{1}/{2}'
                            '/{3}/{4}/{5}/latest/{6}/*'.
                            format(model, experiment, time_frequency, realm,
                                   cmor_table, ensemble, variable)))


def clean_cubelist_atts(cubelist):
    """Wrapper for common_dict - modifies the attributes of a list of
    dictionaries to ensure error-free concatenation."""
    att_dicts = [cube.attributes for cube in cubelist]
    realizations = [d['realization'] for d in att_dicts]
    att_dicts = {k: v for k, v in zip(realizations, att_dicts)}
    for realization in realizations:
        indices = [i for i, j in enumerate(realizations) if j == realization]
        for index in indices:
            cubelist[index].attributes = att_dicts[realization]


def add_realization_number(cube):
    """Add a realization number to a cube's auxiliary coordinates.
    Skip if no realization number is found.
    """
    if 'realization' in cube.attributes.keys():
        realization_number = cube.attributes['realization']
    else:
        realization_number = 1
    realization_number = np.int32(realization_number)
    cube.add_aux_coord(iris.coords.AuxCoord(realization_number,
                       'realization'))

def clean(cube, field, fname):
    guess_bounds(cube)
    add_realization_number(cube)


def guess_bounds(cube):
    """Guess the bounds for the cubes coordinates if there are none."""
    for coord in cube.dim_coords:
        condition = ((coord.bounds is None) and len(coord.points) != 1)
        if condition:
            try:
                coord.guess_bounds()
            except:
                pass


def check_for_timepoint_duplicates(cubes):
    """Remove duplicated time points from a list of cubes so they can
    be concatenated.

    Args:
        cubes (list or iris.cube.CubeList)

    """
    new_cubes = iris.cube.CubeList()
    new_cubes.append(cubes[0])

    # For each cube, append that part of it which doesn't overlap in time
    # with the final cube in the new cube list
    for cube in cubes:
        dates1 = utils.make_datetimes(new_cubes[-1])
        dates2 = utils.make_datetimes(cube)
        non_overlap = np.where(dates2 > dates1[-1])
        if len(non_overlap[0]) != 0:
            new_cubes.append(cube[non_overlap])

    return new_cubes


def check_realizations_timepoint_duplicates(cubes):
    cubes = ensure_time_axis(cubes)
    new_cubes = iris.cube.CubeList()
    # Deal with duplicated data
    realizations = set([cube.coord('realization').points[0]
                        for cube in cubes])
    for realization in realizations:
        constraint = iris.Constraint(coord_values={'realization':
                                                   lambda r: r ==
                                                   realization})
        cubes_realization = cubes.extract(constraint)
        new_cubes.extend(check_for_timepoint_duplicates(cubes_realization))

    return new_cubes


def check_coords(cubes, write_to='./.jazz/offending_cube'):
    """Check coordinates are matching.  If they are not this could be
    quite a problem!  However, some models' have files which read in with
    slightly different coordinates (CCSM4, for example).  In this case
    the difference is miniscule so we can safely replace the coordinates.
    This method replaces coordinates but also informs the user it is doing
    this. It also prints and optionally saves the summary of the offending
    cube.

    Args:
        cubes (iris.cube.CubeList): list of cubes to check
        write_to (Optional[str]): path to which to write warnings

    """
    # Get the names of the spatial coords
    coord_names = [coord.name() for coord in cubes[0].dim_coords]
    if 'time' in coord_names:
        coord_names.remove('time')

    for coord_name in coord_names:
        # Make a list of the coordinates' points for each cube
        points_list = [cube.coord(coord_name).points for cube in cubes]

        # Loop over the list of points for all the cubes
        for p in xrange(len(points_list)-1):

            # If the coordinates are different from the first set,
            # replace them with the first set
            if not (points_list[p+1] == points_list[0]).all():
                cubes[p+1].replace_coord(cubes[0].coord(coord_name))

                # Notify user
                warnings.warn('Replacing the coordinates of a cube. '
                              'Offending cube is {}'.
                              format(cubes[p+1].summary()))

                if write_to is not None:
                    utils.mkdirp(write_to)
                    utils.write_file(cubes[p+1].summary(),
                                     '{0}_{1}_{2}'.format(write_to,
                                                          coord_name,
                                                          p))


def homogenise_air_pressure(cubes):
    """Some files extend to higher altitude.  Extract only those altitude
    points common to all cubes.

    Args:
        cubes (iris.cube.CubeList)

    Returns:
        iris.cube.CubeList

    """
    # Get smallest air pressure coord
    coords = [cube.coord('air_pressure').points for cube in cubes]

    # Find smallest air pressure coord
    smallest_size = np.min([coord.shape[0] for coord in coords])
    smallest_pts = [coord for coord in coords if coord.shape[0] ==
                    smallest_size]

    # Form a constraint based on the smallest points
    constraint = iris.Constraint(coord_values={'air_pressure': lambda p:
                                               p >= np.min(smallest_pts)})
    cubes = cubes.extract(constraint)
    for cube in cubes:
        if len(cube.coord('air_pressure').points) > 1:
            cube.coord('air_pressure').bounds = None
            cube.coord('air_pressure').guess_bounds()

    return cubes


def ensure_time_axis(cubes):
    """Add a time axis if the cube has a singleton (i.e. scalar) time
    dimension.

    Args:
        cubes (iris.cube.CubeList)

    Returns:
        iris.cube.CubeList

    """
    new_list = iris.cube.CubeList()
    for cube in cubes:
        aux_coord_names = [coord.name() for coord in cube.aux_coords]
        if 'time' in aux_coord_names:
            new_list.append(iris.util.new_axis(cube, 'time'))
        else:
            new_list.append(cube)
    return new_list


def fetch(location, constraint=None):
    """Fetch data from a netCDF or a directory containing netCDFs using iris.
    A common cause of concatenate errors is when the time data aren't
    contiguous. This is checked for and corrected.

    Args:
        location (str): path to netCDF file or directory
        constraint (iris.Constraint, optional): constraint for the cubes

    Returns:
        iris.Cube

    """
    # Check if the supplied location is a single file or a directory
    cubes = iris.load(location, constraint, callback=clean)

    # Extract only the data, not the ancillary cubes
    var_name = location.split('/')[-2]
    cubes = cubes.extract(iris.Constraint(cube_func=lambda cube:
                                          cube.var_name == var_name))

    if len(cubes) == 1:
        cube = cubes[0]
    else:
        #clean_cubelist_atts(cubes)
        iris.util.unify_time_units(cubes)
        from iris.experimental.equalise_cubes import equalise_attributes
        equalise_attributes(cubes)
        coord_names = [coord.standard_name for coord in cubes[0].coords()]
        if 'time' in coord_names:
            cubes = check_realizations_timepoint_duplicates(cubes)
        if 'air_pressure' in coord_names:
            cubes = homogenise_air_pressure(cubes)
        check_coords(cubes)

        cubes = cubes.concatenate()
        realizations = [cube.coord('realization').points[0]
                        for cube in cubes]
        if len(set(realizations)) != 1:
            cubes = iris.cube.CubeList(utils.make_common_in_time(*cubes))

        cube = cubes.merge_cube()

    return cube



def available_models(experiments, variables, frequencies, realms,
                     cmor_tables, outfile=None):
    """Get a list of models with the required data.

    Args:
        experiments (list): list of experiments required
        variables (list): list of variables required
        frequencies (list): list of time frequency
        realms (list): CMIP5 realm
        cmor_tables (list): CMOR table
        outfile (Optional[str]): full file path to optionally write out data

    Returns:
        set

    """
    models = []
    iterator = itertools.product(experiments, variables, frequencies,
                                 realms, cmor_tables)
    for experiment, variable, frequency, realm, cmor_table in iterator:
        files = availability('*', experiment, frequency, realm, cmor_table,
                             'r1i1p1', variable)
        models.append([f.split('/')[7] for f in files])

    # Get the common models in the lists of lists
    result = set(models[0])
    for item in models[1:]:
        result.intersection_update(item)

    # Optionally write the file
    if outfile != None:
        utils.write_file(sorted(result), outfile)

    return sorted(result)
