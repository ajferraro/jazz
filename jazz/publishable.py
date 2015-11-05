"""Making publication-ready plots."""

import matplotlib.pyplot as plt
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


plot_styling()
