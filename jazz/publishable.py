"""Making publication-ready plots."""

import matplotlib.pyplot as plt
import os
import string


def plot_styling(mode='jazz_paper'):
    try:
        plt.style.use(mode)
    except:
        print 'No style sheet available'


def label_subplots(ax, case='lower', **kwargs):
    """Put letter labels at top right corner of subplots.

    Args:
        ax (array of matplotlib axes): axes to label

    """
    n_axes = ax.size

    if case is 'upper':
        labels = list(string.ascii_uppercase)[0:n_axes]
    else:
        labels = list(string.ascii_lowercase)[0:n_axes]

    for i, axes in enumerate(ax.flat):
        axes.set_title(labels[i], fontdict=kwargs, loc='left')


def save(path, exts=['png'], close=True, verbose=True):
    """Save a figure from pyplot.
    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.
    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.
    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.

    http://www.jesshamrick.com/2012/09/03/saving-figures-from-pyplot/
    """

    if not isinstance(exts, list):
        exts = list(exts)

    for ext in exts:

        # Extract the directory and filename from the given path
        directory = os.path.split(path)[0]
        filename = "%s.%s" % (os.path.split(path)[1], ext)
        if directory == '':
            directory = '.'

        # If the directory does not exist, create it
        if not os.path.exists(directory):
            os.makedirs(directory)

        # The final path to save to
        savepath = os.path.join(directory, filename)

        if verbose:
            print("Saving figure to '%s'...\n" % savepath),

        # Actually save the figure
        plt.savefig(savepath)

    # Close it
    if close:
        plt.close()

    if verbose:
        print("Done")


def add_colorbar(fig, im, axes, axes2=None, padfrac=0.1, wfrac=0.1, lfrac=0.05,
                 **kw):
    """Add a colour bar based on the dimensions of the given axes.  Unlike
    matplotlib.colorbar.make_axes, this does not 'steal' space from the given
    axes but simply matches their dimensions.

    Args:
        im (mappable): used to generate the colour bar.
        axes (axes): axes to get dimensions from.
        axes2 (Optional[axes]): additional axes if the colorbar should span
            axes.
        padfrac (Optional[float]): fraction of the plot size to pad between
            axes and colour bar.
        wfrac (Optional[float]): width of colour bar in fraction of plot size.
        lfrac (Optional[float]): length of colour bar to subtract in fraction
            of plot size.

    Returns:
        colorbar

    """
    if not 'orientation' in kw:
        kw['orientation'] = 'vertical'

    divider = axes_grid1.make_axes_locatable(axes)
    position = divider.get_position()
    if axes2 is not None:
        divider2 = axes_grid1.make_axes_locatable(axes2)
        position2 = divider.get_position()
        position = (position[0], position[1], position[2]+position2[2], position[3])

    if kw['orientation'] is 'horizontal':
        left = position[0] + lfrac*position[2]
        width = position[2] - 2*lfrac*position[2]
        bottom = position[1] - (padfrac+wfrac)*position[3]
        height = wfrac * position[3]
    elif kw['orientation'] is 'vertical':
        left = position[0] + (1+padfrac)*position[2]
        width = wfrac * position[2]
        bottom = position[1] + lfrac*position[3]
        height = position[3] - lfrac*position[3]

    cax = fig.add_axes((left, bottom, width, height))
    return plt.colorbar(im, cax=cax, **kw)


plot_styling()
