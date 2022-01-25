import numpy as np
import matplotlib.pyplot as plt
import functions as functions
import matplotlib as mlib
viridis_cmap = mlib.cm.get_cmap('viridis')
mlib.rcParams.update({'font.size': 24})


def make_plots(values):
    #remove the first somatic spike
    values -= 1
    values[np.where(values < 0)] = 0
    plt.matshow(values, origin='lower', cmap = plt.cm.Purples)
    cbar = plt.colorbar()
    cbar.set_ticks([])
    plt.axis('off')
    plt.show()


def main():
    soma_frequencies = np.loadtxt("soma_frequency_synapse_big.csv")
    make_plots(soma_frequencies)


if __name__ == '__main__':
    main()
