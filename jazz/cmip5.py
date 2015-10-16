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
    realizations = set([d['realization'] for d in att_dicts])
    att_dicts = [[d for d in att_dicts if d['realization'] == r]
                 for r in realizations]
    att_dicts = [dict_sublist[0] for dict_sublist in att_dicts]
    for cube in cubelist:
        cube.attributes = att_dicts[cube.attributes['realization']-1]


def add_realization_number(cube):
    realization_number = cube.attributes['realization']
    cube.add_aux_coord(iris.coords.AuxCoord(np.int32(realization_number),
                                            'realization'))

def clean(cube, field, fname):
    guess_bounds(cube)
    add_realization_number(cube)


def guess_bounds(cube):
    """Guess the bounds for the cubes coordinates if there are none."""
    for coord in cube.coords():
        if (coord.bounds is None) and (coord.shape != (1,)):
            coord.guess_bounds()


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

    new_cubes = iris.cube.CubeList()
    # Deal with duplicated data
    realizations = set([cube.coord('realization').points[0] for cube in cubes])
    for realization in realizations:
        constraint = iris.Constraint(coord_values={'realization':
                                                   lambda r: r == realization})
        cubes_realization = cubes.extract(constraint)
        new_cubes.extend(check_for_timepoint_duplicates(cubes_realization))

    return new_cubes


def check_coords(cubes, write_to='./offending_cube'):
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

    Returns:
        iris.cube.CubeList

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
                    utils.write_file(cubes[p+1].summary(),
                                     '{0}_{1}_{2}'.format(write_to,
                                                          coord_name,
                                                          p))


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

    if os.path.isfile(location):
        cube = cubes[0]
    else:
        clean_cubelist_atts(cubes)
        iris.util.unify_time_units(cubes) # NEW
        from iris.experimental.equalise_cubes import equalise_attributes
        equalise_attributes(cubes)
        coord_names = [coord.standard_name for coord in cubes[0].coords()]
        if 'time' in coord_names:
            cubes = check_realizations_timepoint_duplicates(cubes)
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





