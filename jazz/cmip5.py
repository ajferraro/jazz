"""
Module for handling CMIP5 files
"""
import glob
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




def common_dict(dict_list):
    """Take a list of dictionaries and delete any keys which have different
    values in the list.
    """
    for compdict in dict_list:

        # Some cubes have the comment attribute, but most don't. Best to 
        # just get rid of it
        if 'comment' in compdict:
            del compdict['comment']

        differences = [k for k,v in dict_list[0].iteritems() if compdict[k]!=v]
    return {k:v for k,v in dict_list[0].iteritems() if k not in differences}


def clean_cubelist_atts(cubelist):
    """Wrapper for common_dict - modifies the attributes of a list of 
    dictionaries to ensure error-free concatenation."""
    attributes = [cube.attributes for cube in cubelist]
    common_attributes = common_dict(attributes)
    for cube in cubelist:
        cube.attributes = common_attributes


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
        if (coord.bounds == None) and (coord.shape != (1,)):
            coord.guess_bounds()


def check_for_timepoint_duplicates(cubes):
    """Remove duplicated time points from a list of cubes so they can 
    be concatenated.

    Args:
        cubes (list or iris.cube.CubeList)

    """
    for i in xrange(len(cubes)-1):
        dates1 = utils.make_datetimes(cubes[i])
        dates2 = utils.make_datetimes(cubes[i+1])
        non_overlap = np.where(dates2 > dates1[-1])
        if len(non_overlap[0]) != 0:
            cubes[i+1] = cubes[i+1][non_overlap]
        else:
            # If all overlapping just remove one of the paired cubes
            cubes.remove(cubes[i+1])
    return cubes


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
        try:
            cubes = cubes.concatenate()
            realizations = [cube.coord('realization').points[0]
                            for cube in cubes]
            if len(set(realizations)) != 1:
                cubes = iris.cube.CubeList(utils.make_common_in_time(*cubes))
        
            cube = cubes.merge_cube()             
        except:
            # If concatenation fails, try the manual approach - forming a new
            # cube.
            print 'Automatic concatenation using iris failed. Trying manual.'
            newdata = []
            newtime = []
            for cube in cubes:
                newdata.extend(cube.data)
                newtime.extend(cube.coord('time').points)

            newdata = np.array(newdata)
            newtime = np.array(newtime)

            # Define the new, concatenated time coordinate
            newtimecoord = iris.coords.DimCoord(newtime,
                                                standard_name=cube.
                                                coord('time').standard_name,
                                                units=cube.coord('time').units)
            
            # Make the new cube
            dim_coords_and_dims = list(enumerate(cube.dim_coords))
            dim_coords_and_dims = [item[::-1] for item in dim_coords_and_dims]
            dim_coords_and_dims[0] = [newtimecoord, 0]
            cube = iris.cube.Cube(newdata,
                                  standard_name=cube.standard_name,
                                  units=cube.units,
                                  attributes=cube.attributes,
                                  dim_coords_and_dims=dim_coords_and_dims)
        
    return cube



def available_models(experiments, variables, frequency, realm, cmor_table,
                     outfile=None):
    """Get a list of models with the required data.

    Args:
        experiments (list): list of experiments required
        variables (list): list of variables required
        frequency (str): time frequency
        realm (str): CMIP5 realm
        cmor_table (str): CMOR table
        outfile (Optional[str]): full file path to optionally write out data
       
    Returns:
        set

    """
    models = []
    for experiment in experiments:
        for variable in variables:
            files = availability('*', experiment, frequency, realm,
                                       cmor_table, 'r1i1p1', variable)
            models.append([f.split('/')[7] for f in files])
    
    # Get the common models in the lists of lists
    result = set(models[0])
    for item in models[1:]:
        result.intersection_update(item)

    # Optionally write the file
    if outfile != None:
        utils.write_file(sorted(result), outfile)
    
    return sorted(result)





