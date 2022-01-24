#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlib
from matplotlib.lines import Line2D

# define colourmap
viridis_cmap = mlib.cm.get_cmap('viridis')
blues_cmap = mlib.cm.get_cmap('Greens')
purd_cmap = mlib.cm.get_cmap('PuRd')
mlib.rcParams.update({'font.size': 24})

def figure_A():
    plot_data = np.loadtxt("figure4A_data.txt")
    inputs = [32, 35, 38, 42, 45]
    plots_colours = np.linspace(0.2, 0.8, 5)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    count = 0
    for input in inputs:
        ax1.plot(plot_data[count,:], color=purd_cmap(0.99), alpha=plots_colours[count])
        count += 1

    ax1.set_yticks([-40, 20])
    ax1.set_xticks([0, 500])
    ax1.set_xticklabels([0, 5])
    custom_lines = [Line2D([0], [0], color=purd_cmap(0.99), lw=2, alpha = plots_colours[0]),
                    Line2D([0], [0], color=purd_cmap(0.99), lw=2, alpha = plots_colours[1]),
                    Line2D([0], [0], color=purd_cmap(0.99), lw=2, alpha = plots_colours[2]),
                    Line2D([0], [0], color=purd_cmap(0.99), lw=2, alpha = plots_colours[3]),
                    Line2D([0], [0], color=purd_cmap(0.99), lw=2, alpha = plots_colours[4])]
    ax1.legend(custom_lines, ['32 nA', '35 nA', '38 nA', '41 nA', '44 nA'], frameon=False)
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Membrane voltage (mV)')


def figure_B():
    line_data = np.loadtxt("figure4B_lineData.txt")
    scatter_data = np.loadtxt("figure4B_scatterData.txt")

    inputs_scatter = [32, 35, 38, 42, 45]
    inputs = np.linspace(25, 50, 50)
    plots_colours = np.linspace(0.2, 0.8, 5)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.plot(line_data[0,:], line_data[1,:], color = purd_cmap(0.99), alpha = 0.15, linewidth = 3.5)
    for count in range(0, 5):
        ax1.scatter(scatter_data[0, count], scatter_data[1, count], c = purd_cmap(0.99), alpha = plots_colours[count])
    ax1.set_yticks([-40, 20])
    ax1.set_xticks([25, 50])
    ax1.set_xlabel('Injected current (nA)')
    ax1.set_ylabel('dCaAP amplitude (mV)')

    
def figure_C():
    plot_data = np.loadtxt("figure4C_data.txt")
    inputs = [32, 35, 38, 42, 45]
    plots_colours = np.linspace(0.2, 0.8, 5)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    count = 0
    for input in inputs:
        ax1.plot(plot_data[count, :], color=viridis_cmap(0.15), alpha = plots_colours[count])
        count += 1

    custom_lines = [Line2D([0], [0], color=viridis_cmap(0.15), lw=2, alpha = plots_colours[0]),
                    Line2D([0], [0], color=viridis_cmap(0.15), lw=2, alpha = plots_colours[1]),
                    Line2D([0], [0], color=viridis_cmap(0.15), lw=2, alpha = plots_colours[2]),
                    Line2D([0], [0], color=viridis_cmap(0.15), lw=2, alpha = plots_colours[3]),
                    Line2D([0], [0], color=viridis_cmap(0.15), lw=2, alpha = plots_colours[4])]
    # mlib.rcParams.update({'font.size': 13})
    ax1.legend(custom_lines, ['32 nA', '35 nA', '38 nA', '41 nA', '44 nA'], frameon=False)
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Membrane voltage (mV)')
    ax1.set_yticks([-70, 20])
    ax1.set_xticks([0, 500])
    ax1.set_xticklabels([0, 5])
    

def figure_D():
    line_data = np.loadtxt("figure4D_lineData.txt")
    scatter_data = np.loadtxt("figure4D_scatterData.txt")
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    plots_colours = np.linspace(0.2, 0.8, 5)
    ax1.plot(line_data[0,:], line_data[1,:], color=viridis_cmap(0.15), alpha=0.15, linewidth=3.5)
    for count in range(0, 5):
        ax1.scatter(scatter_data[0, count], scatter_data[1, count], c=viridis_cmap(0.15), alpha=plots_colours[count])
    ax1.set_yticks([-60, 20])
    ax1.set_xticks([25, 50])
    ax1.set_xlabel('Injected current (nA)')
    ax1.set_ylabel('somatic AP amplitude (mV)')

def main():
    figure_A()
    figure_B()
    figure_C()
    figure_D()
    plt.show()


if __name__ == '__main__':
    main()
