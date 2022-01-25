import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlib
from matplotlib.lines import Line2D

# define colourmap
viridis_cmap = mlib.cm.get_cmap('viridis')
blues_cmap = mlib.cm.get_cmap('Greens')
purd_cmap = mlib.cm.get_cmap('PuRd')
mlib.rcParams.update({'font.size': 24})

def figure_AB():
    plot_data = np.loadtxt("dataFiles/figure3AB_data.txt")
    timesteps = np.linspace(0, int(1000), int(10000))

    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')

    ax1.plot(timesteps[2000:], plot_data[1,:], color=purd_cmap(0.85), alpha=1, label='$V_{dend}$', linewidth=2.5)
    ax1.plot(timesteps[2000:], plot_data[0,:], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth = 2.5)
    ax1.legend(frameon=False)
    ax1.set_yticks([-75, 20])
    ax1.set_xticks([200, 1000])
    ax1.set_xticklabels(['0', '800'])
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Membrane voltage (mV)')


    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')

    ax1.plot(timesteps[2000:2500], plot_data[0, 0:500], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth = 2.5)
    ax1.plot(timesteps[2000:2500], plot_data[1, 0:500], color=purd_cmap(0.8), alpha=1, label='$V_{dend}$', linewidth = 2.5)
    ax1.legend(frameon=False)
    ax1.set_ylabel("Membrane voltage (mV)")
    ax1.set_xlabel('Time (ms)')
    ax1.set_yticks([-80, 20])
    ax1.set_xticks([200, 250])
    ax1.set_xticklabels(['0', '50'])


def figure_C():
    plot_data = np.loadtxt("dataFiles/figure3C_data.txt")
    fig = plt.figure()
    total_time = 3
    time = total_time * 1000
    timesteps = np.linspace(0, int(time), int(time*10))

    for count in [0,1,2]:
        ax1 = plt.subplot(3,1, count+1)
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.plot(timesteps[:], plot_data[count, :], color=purd_cmap(0.85), alpha=1, label='$V_{dend}$', linewidth=1.5)
        ax1.set_yticks([-75, 20])
        ax1.set_xticks([0, 3000])
        if count == 0:
            ax1.legend(frameon=False, loc=1)
        count += 1

    ax1.set_ylabel('Membrane voltage (mV)')
    ax1.set_xlabel('Time (ms)')


def figure_DE():
    plot_data = np.loadtxt("dataFiles/figure3DE_data.txt")
    timesteps = np.linspace(0, int(0.1*1000), int(0.1*10000))
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')

    ax1.plot(timesteps[150:800], plot_data[1,:], color=purd_cmap(0.8), alpha=1, label='$V_{dend}$', linewidth=2.5)
    ax1.plot(timesteps[150:800], plot_data[0,:], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth = 2.5)
    ax1.legend(frameon=False, loc=1)
    ax1.set_yticks([-75, 20])
    ax1.set_xticks([15, 80])
    ax1.set_xticklabels(['0', '80'])
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Membrane voltage (mV)')


    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.plot(timesteps[150:250], plot_data[0, :100], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth = 2.5)
    ax1.plot(timesteps[150:250], plot_data[1, :100], color=purd_cmap(0.8), alpha=1, label='$V_{dend}$', linewidth = 2.5)
    ax1.legend(frameon=False)
    ax1.set_ylabel("Membrane voltage (mV)")
    ax1.set_xlabel('Time (ms)')
    ax1.set_yticks([-75, 20])
    ax1.set_xticks([15, 25])
    ax1.set_xticklabels(['0', '10'])


def figure_F():
    plot_data = np.loadtxt("dataFiles/figure3F_data.txt")
    fig = plt.figure()
    total_time = 0.3
    time = total_time * 1000
    timesteps = np.linspace(0, int(time), int(time*10))

    for count in [0,1,2]:
        ax1 = plt.subplot(3,1, count+1)
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.plot(timesteps[:], plot_data[count, :], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth = 1.5)
        ax1.set_yticks([-75, 20])
        ax1.set_xticks([0, 300])
        if count == 0:
            ax1.legend(frameon=False, loc=1)
        count += 1

    ax1.set_ylabel('Membrane voltage (mV)')
    ax1.set_xlabel('Time (ms)')


def main():
    figure_AB()
    figure_C()
    figure_DE()
    figure_F()
    plt.show()


if __name__ == '__main__':
    main()
