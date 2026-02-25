import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from matplotlib.lines import Line2D

# da/dt = a + b - 5
#
# db/dt = 2a - 8b

da = 5
db = 5

a_fp = 4
b_fp = 1

def phase_portrait(system,
                    x0=None,
                    xlim=None,
                    ylim=None,
                    dx=5,
                    dy=5,
                    grid_points=100,
                    figsize=(6, 6),
                    fontsize=12,
                    streamplot_kwargs=None,
                    plot_nullclines=True,
                    title=None,
                    xlabel=None,
                    ylabel=None,
                    nullcline_labels=None,
                    fixedpoint_label='fixed point',
                    show_legend=True):
    """Create a phase portrait for a 2D system.

    Parameters
    ----------
    system : callable
        Function that accepts a length-2 array-like [x, y] (or two arrays) and
        returns [dx/dt, dy/dt]. The function should support numpy arrays.
    x0 : array-like, optional
        Initial guess for finding a fixed point with fsolve. Default [0, 0].
    xlim, ylim : tuple, optional
        Axis limits as (min, max). If None they are set relative to the fixed
        point using dx/dy.
    dx, dy : float
        Range around the fixed point to build the plotting window when xlim/
        ylim are not provided.
    grid_points : int
        Number of points along each axis for the vector field grid.
    figsize : tuple
        Figure size passed to plt.subplots.
    fontsize : int
        Font size for labels and title.
    streamplot_kwargs : dict, optional
        Additional keyword args forwarded to ax.streamplot.
    plot_nullclines : bool
        If True, draw contour lines where dx/dt=0 and dy/dt=0.

    Returns
    -------
    fig, ax
    """
    # Prepare default kwargs
    if streamplot_kwargs is None:
        streamplot_kwargs = {'density': 1.0, 'color': 'gray'}

    # Find fixed point
    x0 = np.asarray(x0) if x0 is not None else np.array([0.0, 0.0])
    try:
        fp = fsolve(system, x0)
    except Exception:
        # fallback to zero if fsolve fails
        fp = np.array([0.0, 0.0])

    # Determine plotting ranges
    if xlim is None:
        x_min, x_max = fp[0] - dx, fp[0] + dx
    else:
        x_min, x_max = xlim
    if ylim is None:
        y_min, y_max = fp[1] - dy, fp[1] + dy
    else:
        y_min, y_max = ylim

    x = np.linspace(x_min, x_max, int(grid_points))
    y = np.linspace(y_min, y_max, int(grid_points))
    X, Y = np.meshgrid(x, y)

    # Evaluate vector field. Accept systems that return lists or arrays.
    d = system([X, Y])
    # ensure we have two arrays
    dX = np.asarray(d[0])
    dY = np.asarray(d[1])

    fig, ax = plt.subplots(figsize=figsize)

    # Streamplot
    ax.streamplot(X, Y, dX, dY, **streamplot_kwargs)

    # Plot fixed point
    ax.plot(fp[0], fp[1], marker='x', color='k', markersize=8, markeredgewidth=2)
    if fixedpoint_label:
        ax.annotate(fixedpoint_label, xy=(fp[0], fp[1]), xytext=(5, 5), textcoords='offset points', fontsize=fontsize)

    # Nullclines via contour at level 0
    handles = []
    labels = []
    if plot_nullclines:
        try:
            c1 = ax.contour(X, Y, dX, levels=[0], colors='C0', linestyles='--', linewidths=1)
            c2 = ax.contour(X, Y, dY, levels=[0], colors='C1', linestyles='-.', linewidths=1)
            # create legend handles
            # use provided labels if any
            if nullcline_labels is None:
                nc_labels = ('dx/dt = 0', 'dy/dt = 0')
            else:
                # ensure we have at least two labels
                nc_labels = tuple(nullcline_labels)
                if len(nc_labels) < 2:
                    defaults = ('dx/dt = 0', 'dy/dt = 0')
                    nc_labels = tuple(list(nc_labels) + list(defaults[len(nc_labels):]))
            handles.append(Line2D([0], [0], color='C0', linestyle='--'))
            labels.append(nc_labels[0])
            handles.append(Line2D([0], [0], color='C1', linestyle='-.'))
            labels.append(nc_labels[1])
        except Exception:
            # if contouring fails (e.g., non-finite values), skip nullclines
            pass

    # axis labels and formatting
    ax.set_xlabel(xlabel if xlabel is not None else 'x', fontsize=fontsize)
    ax.set_ylabel(ylabel if ylabel is not None else 'y', fontsize=fontsize)
    if title is not None:
        ax.set_title(title, fontsize=fontsize)
    else:
        ax.set_title('Phase portrait', fontsize=fontsize)

    # Add a legend combining nullcline proxies (if any)
    if show_legend and handles:
        # include fixed point in legend too if requested
        fp_label = fixedpoint_label if fixedpoint_label else None
        if fp_label:
            handles = [Line2D([0], [0], color='k', marker='x', linestyle='')] + handles
            labels = [fp_label] + labels
        ax.legend(handles, labels, fontsize=fontsize)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    return fig, ax

alpha, beta, gamma, delta = 2, 1.1, 1, 0.9
a0, b0 = 1, 0.5
def dadt(a, b):
    return alpha*a - beta*a*b #a**2 + b - b**2 - 5

def dbdt(a, b):
    return -gamma*b + delta*a*b #a**2 - b**2 - 8

def nca(alpha, beta):
    return alpha/beta

def ncb(gamma, delta):
    return gamma/delta

def system(vars):
    a, b = vars[0], vars[1]
    return [dadt(a, b), dbdt(a, b)]




# Example usage: create a phase portrait for the example linear system
if __name__ == '__main__':
    fig, ax = phase_portrait(system,
                             x0=[a0, b0],
                             dx=5,
                             dy=5,
                             grid_points=200,
                             figsize=(6, 6),
                             fontsize=12,
                             streamplot_kwargs={'density': 1.5, 'color': 'gray'})
    plt.show()
