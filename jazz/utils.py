import numpy as np
from scipy import stats
import datetime

import iris

import cmip5


def anomalies(cube, frac=False, clim=None):
    """Calculate anomalies from a monthly climatology.  Can handle non-January
    start months and cubes of varying dimensions, as long as the time
    dimension is first.

    Args:
        cube (iris.Cube): cube from which to calculate anomalies.
        frac (bool, optional): whether to calculate anomalies as %
            deviations from the climatology
        climatology: placeholder for the climatology data - if this is passed
            as a numpy array it can be used as an extra optional output

    Returns
        iris.Cube

    """
    anom_cube = cube.copy()
    ntim = cube.shape[0]
    n_months_in_yr = 12
    non_time_dims = cube.data.shape[1:]
    climatology = np.zeros((n_months_in_yr,)+non_time_dims)

    # Calculate dates for the cube
    dates = cube.coord('time').units.num2date(cube.coord('time').points)
    months = np.array([d.month for d in dates])
    startmonth = months[0]
    
    # Make the climatology
    for m in xrange(n_months_in_yr):
        month_indices = np.where(months == m+1)
        climatology[m] = cube[month_indices].collapsed('time', 
                                                       iris.analysis.MEAN).data

    # Calculate anomalies
    for t in xrange(ntim):
        m = (t+startmonth)%12
        anom_cube.data[t] = cube.data[t] - climatology[m-1]

    if frac == True:
        anom_cube = 100*(anom_cube/cube.collapsed('time', iris.analysis.MEAN))
    
    if clim is not None:
        clim[:] = climatology[0:12]

    return anom_cube


def detrended_anomalies(cube):
    """Shortcut function to calculate anomalies of a detrended cube."""
    return anomalies(detrend(cube))


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


def time_mean(*args):
    """Calculate the time average of iris CUbes"""
    ret = []
    for cube in args:
        ret.append(cube.collapsed('time', iris.analysis.MEAN))

    if len(ret) == 1:
        ret = ret[0]
    return ret


def area_weighted(cube):
    """Weight a cube by latitude."""
    grid_areas = iris.analysis.cartography.area_weights(cube)
    return cube * grid_areas / np.mean(grid_areas)


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


# Linear horizontal regrid if the supplied regridding coordinates are 
# different to the cube's original coordinates
def regrid(cube, gridpts):
    if not (np.array_equal(gridpts[0][1], cube.coord('latitude').points) & 
            np.array_equal(gridpts[1][1], cube.coord('longitude').points)):
        cube = iris.analysis.interpolate.linear(cube,gridpts)
        cmip5.guess_bounds(cube)
    return cube


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


def loop_dictionary(dictionary, function, *args):
    for key,val in dictionary.items():
        dictionary[key] = function(val, *args)
    return dictionary


def detrend(cube):
    """Linearly detrend data in an iris cube.  The time dimension must be 
    first.

    Args:
        cube (iris.Cube)

    Returns:
        iris.Cube

    """
    shape = cube.shape
    # Reshape so the time dimension is followed by a single non-time dimension
    data = cube.data.reshape(shape[0], np.prod(shape[1:]))
    newshape = data.shape
    time = np.arange(shape[0])
    trendline = np.zeros_like(data)
    for i in xrange(newshape[1]):
        regresults = stats.mstats.linregress(time, data[:,i])
        slope = regresults[0]
        intercept = regresults[1]
        trendline[:,i] = slope*time
        
    newdata = (data-trendline).reshape(shape)
    return cube.copy(data=newdata)


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



def make_datetimes(cube):
    pseudo = cube.coord('time').units.num2date(cube.coord('time').points)
    return np.array([datetime.datetime(d.year, d.month, d.day) 
                     for d in pseudo])


def read_file(filename):
    contents = []
    with open(filename, 'r') as file:
        for line in file:
            if not line.strip().startswith("#"):
                contents.append(line.strip('\n'))
    return contents


def write_file(data, filename):
    with open(filename, 'w') as file:
        for line in data:
            file.write('{}\n'.format(line))


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
