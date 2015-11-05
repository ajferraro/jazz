import numpy as np
import os
from scipy import stats
import datetime

import iris
import iris.coord_categorisation as cat

import cmip5


def anomalies(cube, kind='month'):
    """Calculate anomalies from a monthly or yearly climatology.

    Args:
        cube (iris.cube.Cube)
        kind (Optional[str]): 'month' or 'year'

    Returns:
        iris.cube.Cube

    """
    # Calculate climatology
    clim = climatology(cube, kind=kind)

    # Repeat the climatology in time so it can be subtracted from the cube
    ntim_clim = clim.coord('time').points.shape[0]
    ntim_cube = cube.coord('time').points.shape[0]
    nrepeats = ntim_cube/ntim_clim
    coord_name_tuple = [coord.name() for coord in clim.dim_coords]

    # Find out where the time axis is
    time_axis = coord_name_tuple.index('time')
    clim_data = np.expand_dims(clim.data, time_axis)
    clim_data = np.repeat(clim_data, nrepeats,
                          axis=time_axis).reshape(cube.shape)

    return cube.copy(data=cube.data-clim_data)


def annual_mean(cube):
    """Calculate the annual mean from a monthly cube.  The time dimension
    of the new cube takes the month of July to identify each year.

    Args:
        cube (iris.cube.Cube)

    Returns:
        newcube (iris.cube.Cube)

    """
    dates = cube.coord('time').units.num2date(cube.coord('time').points)
    months = np.array([d.month for d in dates])
    where_jans = np.where(months == 1)[0]
    where_juns = np.where(months == 6)[0]

    # Create a new time coordinate, including bounds
    new_time_coord = cube.coord('time')[where_juns]
    bounds = np.zeros(new_time_coord.shape+(2,))
    for y in xrange(bounds.shape[0]):
        bounds[y,0] = y*365 + 1
        bounds[y,1] = y*365 + 365

    new_time_coord.bounds = bounds

    # Calculate annual mean data
    data = np.array([cube[i:i+12].collapsed('time', iris.analysis.MEAN).data
                     for i in where_jans])

    # Make the new cube
    dim_coords_and_dims = list(enumerate(cube.dim_coords))
    dim_coords_and_dims = [item[::-1] for item in dim_coords_and_dims]
    dim_coords_and_dims[0] = (new_time_coord, 0)
    newcube = iris.cube.Cube(data, standard_name=cube.standard_name,
                             long_name=cube.long_name, units=cube.units,
                             dim_coords_and_dims=dim_coords_and_dims)
    return newcube


def area_mean(*args):
    """Calculate the area-weighted spatial mean of iris Cubes."""
    ret = []
    for cube in args:
        grid_areas = iris.analysis.cartography.area_weights(cube)
        ret.append(cube.collapsed(['latitude', 'longitude'],
                                  iris.analysis.MEAN, weights=grid_areas))
    if len(ret) == 1:
        ret = ret[0]
    return ret


def area_weighted(cube):
    """Weight a cube by latitude."""
    grid_areas = iris.analysis.cartography.area_weights(cube)
    return cube * grid_areas / np.mean(grid_areas)


def climatology(cube, kind='month'):
    """Calculate a climatology for a cube.  Can do monthly or yearly.

    Args:
        cube (iris.cube.Cube)
        kind (Optional[str]): 'month' or 'year'

    Returns:
        iris.cube.Cube

    """
    aux_coords = [aux_coord.name() for aux_coord in cube.aux_coords]
    if 'year' not in aux_coords:
        cat.add_year(cube, 'time')
    if 'month' not in aux_coords:
        cat.add_month(cube, 'time')
    return cube.aggregated_by(kind, iris.analysis.MEAN)


def common_grid(ref_index, *args):
    """Put all supplied cubes on a common horizontal grid.

    Args:
        ref_index (int): index of the cube which is to be used as the
            reference cube, i.e. the one which supplies the common grid
        *args (iris.Cube): iris cubes

    Returns:
         iris.Cube [same number as supplied by *args]

    """
    if len(args) < 2:
        raise Exception, ('Need to supply minimum of two cubes so they can'
                          ' be put on a common grid')

    # Name the reference cube
    ref_cube = args[ref_index]

    # Get the horizontal grid of the reference cube
    gridpts = [('latitude', ref_cube.coord('latitude').points),
               ('longitude', ref_cube.coord('longitude').points)]

    # Interpolate all supplied cubes to the reference horizontal grid
    cubes = [regrid(cube, gridpts) for cube in args]

    return cubes


def constrain_on_time(daterange, *args):
    """Extract the part of the cube falling in the time interval given by
    `daterange'.

    Args:
        daterange (2-element list of datetime.datetime)
        *args: any number of (iris.cube.Cube)

    Returns:
        any number of iris.cube.Cube

    """
    # Check for midnight in daterange. This causes problems since it is
    # technically in two different days (or none!). Add an hour on to
    # account for this.
    daterange = [d if d.hour != 0 else d + datetime.timedelta(seconds=3600)
                 for d in daterange]

    dates = [args[0].coord('time').units.date2num(d) for d in daterange]
    constraint = iris.Constraint(time= lambda c: dates[0] <= c <= dates[1])

    out = [cube.extract(constraint) for cube in args]
    if len(out) == 1:
        out = out[0]

    return out


def constrain_dates(daterange, *data):
    """Constrain cubes to cover a specified daterange.

    Args:
        daterange (str): time range covered (can be generated from
            utils.gen_timestring())
        *args (iris.Cube): cubes to constrain

    Returns:
        list of iris.Cube, or if a single cube is provided in *data, a single
            iris.Cube.

    """
    # Get datetime objects representing the daterange
    d1, d2 = [datetime.datetime(int(d[0:4]), int(d[4:6]), 15)
              for d in (daterange.split('-'))]

    ret = []
    for cube in data:
        # Get the datetime object array for each cube's time dimension
        cubecal = cube.coord('time').units.num2date(cube.coord('time').points)
        cubecal = np.array([datetime.datetime(d.year, d.month, 15)
                            for d in cubecal])
        # Find indices of datetime array within date range
        mask = (cubecal >= d1) & (cubecal <= d2)
        # Add masked cube to output array
        ret.append(cube[np.where(mask==True)[0]])

    # If only a single cube was passed in *data we don't need to return a list
    if len(ret) == 0:
        ret = ret[0]
    return ret


def detrend(cube):
    """Linearly detrend data in an iris cube.  Uses stats.mstats which allows
    for masked data points.

    Args:
        cube (iris.Cube)

    Returns:
        iris.Cube

    """
    coord_name_tuple = [coord.name() for coord in cube.dim_coords]
    time_axis = coord_name_tuple.index('time')
    new_data = detrend_array(cube.data, axis=time_axis)
    return cube.copy(data=new_data)


def detrended_anomalies(cube):
    """Shortcut function to calculate anomalies of a detrended cube."""
    return anomalies(detrend(cube))


def detrend_array(data, axis=0):
    """This is substantially based on scipy.signal.detrend, but it differs
    in two ways: 1) it removes just the slope rather than the slope and
    the intercept, and 2) it uses stats.mstats.linregress and consequently
    can handle masked data, though this is slower than the matrix-based
    approach of scipy.signal.detrend.

    Args:
        data (array_like): the input data
        axis (Optional[int]): axis along which to detrend

    Returns:
        array_like

    """
    dshape = data.shape

    # Restructure data so that axis is along first dimension and
    #  all other dimensions are collapsed into second dimension.
    rank = len(dshape)
    ntim = dshape[axis]
    newdims = np.r_[axis, 0:axis, axis+1:rank]
    newdata = np.reshape(np.transpose(data, tuple(newdims)),
                         (ntim, np.prod(dshape, axis=0) // ntim))
    newdata = newdata.copy()
    time = np.arange(ntim)
    trend_arr = np.zeros_like(newdata)
    N = newdata.shape[1]

    # Calculate trend and remove it.
    for i in xrange(newdata.shape[1]):
        print (float(i)/N)*100
        regresults = stats.mstats.linregress(time, newdata[:, i])
        trend_arr[:, i] = time*regresults[0]

    newdata = newdata - trend_arr

    # Put data back in original shape.
    tdshape = np.take(dshape, newdims, 0)
    ret = np.reshape(newdata, tuple(tdshape))
    vals = list(range(1, rank))
    olddims = vals[:axis] + [0] + vals[axis:]
    ret = np.transpose(ret, tuple(olddims))
    return ret


def gen_timestring(cube):
    """Generate a string describing time range covered by a cube.

    Args:
        cube (iris.Cube)

    Returns:
        str

    """
    calendar = cube.coord('time').units.num2date(cube.coord('time').points)
    date_start = calendar[0]
    date_end = calendar[-1]
    return '-'.join([str(date_start.year)+str(date_start.month).zfill(2),
                     str(date_end.year)+str(date_end.month).zfill(2)])


def loop_dictionary(dictionary, function, *args):
    """Apply a function to all values in a dictionary.

    Args:
        dictionary (dict): input dictionary
        function (function): function to apply
        *args: arguments to `function`

    Returns:
        dict

    """
    for key,val in dictionary.items():
        dictionary[key] = function(val, *args)
    return dictionary


def make_common_in_time(*args):
    """Subet all cubes so they have the same time dimensions.

    Args:
        *args (iris.Cube)

    Returns:
        list

    """
    # Get the calendar arrays for all cubes
    calendars = [cube.coord('time').units.num2date(cube.coord('time').points)
                 for cube in args]
    # Make all days 16 - we don't care for monthly data but we need to
    # make it consistent for the subsequent comparisons
    calendars = [[datetime.datetime(d.year, d.month, 16) for d in calendar]
                 for calendar in calendars]

    latest_start = max([calendar[0] for calendar in calendars])
    earliest_end = min([calendar[-1] for calendar in calendars])
    out = []
    for cube, calendar in zip(args, calendars):
        i_start = np.where(np.array(calendar) == latest_start)[0]
        i_end = np.where(np.array(calendar) == earliest_end)[0]
        out.append(cube[i_start:i_end+1])

    return out


def make_datetimes(cube):
    """Return an array of datetime.datetime objects correponsing to a cube's
    time coordinate.

    Args:
        cube (iris.cube.Cube): input data

    Returns:
        array of datetime.datetime

    """
    pseudo = cube.coord('time').units.num2date(cube.coord('time').points)
    return np.array([datetime.datetime(d.year, d.month, d.day)
                     for d in pseudo])


def read_file(filename):
    """Read in a text file.

    Args:
        filename (str): file to read

    Returns:
        list

    """
    contents = []
    with open(filename, 'r') as file:
        for line in file:
            if not line.strip().startswith("#"):
                contents.append(line.strip('\n'))
    return contents


def regrid(cube, gridpts):
    """Linear horizontal regrid if the supplied regridding coordinates are
    different to the cube's original coordinates.

    Args:
        cube (iris.cube.Cube): cube to regrid
        gridpts: grid to regrid to

    Returns:
        iris.cube.Cube

    """
    if not (np.array_equal(gridpts[0][1], cube.coord('latitude').points) &
            np.array_equal(gridpts[1][1], cube.coord('longitude').points)):
        cube = iris.analysis.interpolate.linear(cube,gridpts)
        cmip5.guess_bounds(cube)
    return cube


def regrid_data(grid_cube, *regrid_cubes):
    """Regrid arbitrary number of cubes.  Masked data is handled by applying
    a mask if more than half the data points within the target cell are masked.

    Args:
        grid_cube (iris.Cube): reference cube supplying new grid
        *regrid_cubes (iris.Cube): cubes to regrid

    Returns:
        iris.Cube if only one cube regridded, list otherwise

    """
    import iris.experimental.regrid as regrid
    ret_data = []
    rg_func = regrid.regrid_area_weighted_rectilinear_src_and_grid
    for cube in regrid_cubes:
        ret_data.append(rg_func(cube, grid_cube, mdtol=0.5))
    if len(regrid_cubes) == 1:
        ret_data = ret_data[0]
    return ret_data


def spatial_constraint(lat=[-90, 90], lon=[0, 360], inc_equal=True):
    """Make a cube constraint in latitude and/or longitude.

    Args:
        lat (Optional[list]): latitude bounds
        lon (Optional[list]): longitude bounds
        inc_equal (bool): include those points where the latitude or
            longitude is equal to the bounds

    Returns:
        iris.Constraint

    """
    if inc_equal:
        lat_func = lambda l: lat[0] <= l <= lat[1]
        lon_func = lambda l: lon[0] <= l <= lon[1]
    else:
        lat_func = lambda l: lat[0] < l < lat[1]
        lon_func = lambda l: lon[0] < l < lon[1]

    return iris.Constraint(coord_values={'latitude': lat_func,
                                         'longitude': lon_func})


def time_mean(*args):
    """Calculate the time average of iris CUbes"""
    ret = []
    for cube in args:
        ret.append(cube.collapsed('time', iris.analysis.MEAN))

    if len(ret) == 1:
        ret = ret[0]
    return ret


def write_file(data, filename):
    """Write data to a text file."""
    with open(filename, 'w') as file:
        for line in data:
            file.write('{}\n'.format(line))


def mkdirp(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
