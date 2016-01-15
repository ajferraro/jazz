"""Making publication-ready plots."""

import matplotlib.pyplot as plt
import os
import string


def plot_styling(mode='jazz_paper'):
    try:
        plt.style.use(mode)
    except:
        print 'No style sheet available'


def label_subplots(ax, case='lower', xpos=0, ypos=1.05, **kwargs):
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
        axes.text(xpos, ypos, labels[i], fontdict=kwargs,
                  transform=axes.transAxes)


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
            print("Saving figure to '%s'..." % savepath),

        # Actually save the figure
        plt.savefig(savepath)

    # Close it
    if close:
        plt.close()

    if verbose:
        print("Done")

plot_styling()
