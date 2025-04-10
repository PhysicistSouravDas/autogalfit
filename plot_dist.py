"""
.. module:: plot_dist

:Synopsis: Plots the histograms of various parameters and saved.
:Author: Sourav Das
:Year: 2023
"""

from matplotlib import pyplot as plt

def plot_dist(x: list, color, ec, title, xlabel, ylabel, savedir, name, bins=20, density=True):
    """
    Plots a histogram of array x, which contains integers or floats

    x(list): Array of data for which histogram to be plotted
    color(str): Color of the bars in histogram
    ec(str): Edge color of the bars in histogram
    title(str): Title of the plot
    xlabel(str): Label for x-axis
    ylabel(str): Label for y-axis
    savedir(str): Directory for saving the plot
    name(str): Name of the plot to be saved (e.g: Name.png)
    bins(int): No. of bins of the histogram (Default 20)
    density(bool): Whether to normalize the plot (Default True)

    Saves the plots in corresponding directory with plot name and returns
    the status of saving plot.
    """
    #4287f5
    plt.title(f"{title}")
    plt.ylabel(f"{ylabel}")
    plt.xlabel(f"{xlabel}")
    plt.hist(x, bins=bins, color=color, ec=ec, density=density)
    plt.savefig(f'{savedir}/{name}.png')
    plt.clf()

    return f"Plot saved to {savedir} as {name}.png"
