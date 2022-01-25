import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mlib

viridis_cmap = mlib.cm.get_cmap('viridis')
mlib.rcParams.update({'font.size': 24})

def comp_neuron_with_spin():
    # load data from NEURON
    with open('dataFiles/NEURON_t_data.txt') as f:
        neuronTimeArray = f.read().splitlines()
        for count in range(0, len(neuronTimeArray)):
            neuronTimeArray[count] = float(neuronTimeArray[count])
    with open('dataFiles/NEURON_V_data.txt') as f:
        neuronVArray = f.read().splitlines()
        for count in range(0, len(neuronVArray)):
            neuronVArray[count] = float(neuronVArray[count])

    # load data from SpiNNaker
    with open("dataFiles/SpiNNaker_data.txt") as file:
        spinnaker = file.read().splitlines()
        spinnaker = list(filter(None, spinnaker))
        spinnaker = [float(i) for i in spinnaker]

    # load data from Jib2
    jib2_array = []
    with open("dataFiles/Jib2_data.txt") as file:
        voltageArray = file.read().splitlines()
        for i, str_v, in enumerate(voltageArray):
            replacing = 'exit_res_2[{}] = '.format(i)
            if str_v == "":
                break
            jib2_array.append(float(str_v.replace(replacing, '')))


    # calculate errors
    error_jib2=[]
    error_spiNNaker = []
    count=0
    for voltage in neuronVArray[:10001]:
        error_jib2.append(voltage - jib2_array[count])
        error_spiNNaker.append(voltage - spinnaker[count])
        # the lines below compare spike times between NEURON, jib2 and SpiNNaker
        # if neuronVArray[count] > -20 and neuronVArray[count-1] < -20:
        #     print(count)
        # if jib2_array[count] > -20 and jib2_array[count-1] < -20:
        #     print(count)
        # if spinnaker[count] > -20 and spinnaker[count-1] < -20:
        #     print(count)
        count+=1

    # here we print the maximum error values between Jib2&NEURON and SpiNNaker&NEURON
    # print(max(error_spiNNaker))
    # print(max(error_jib2))

    # plot the results
    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    dt = 0.1
    time = 600
    t = time / dt

    plotT = np.linspace(0, time, int(t))

    startT = 0
    endT = 60
    ax1.plot(plotT[startT:endT], neuronVArray[startT:endT], alpha=.7, label='$NEURON$', color=viridis_cmap(1))
    ax1.plot(plotT[startT:endT], spinnaker[startT:endT], alpha=.7, label='$SpiNNaker$', color=viridis_cmap(0.6))
    ax1.plot(plotT[startT:endT], jib2_array[startT:endT], alpha=.7, label='$Jib2$', color=viridis_cmap(0.2))
    ax1.set_ylabel('Membrane voltage (mV)')
    ax1.legend(frameon=False)
    ax1.set_yticks([-75, 25])
    ax1.set_xticks([0, 6])

    ax1 = plt.subplot(223)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.plot(plotT[startT:endT], error_jib2[startT:endT], color=viridis_cmap(0.5), alpha=0.7, label = '$Jib2&NEURON$')
    ax1.plot(plotT[startT:endT], error_spiNNaker[startT:endT], color=viridis_cmap(1), alpha=0.7, label = '$SpiNNaker&NEURON$')
    ax1.legend(frameon=False)
    ax1.set_xticks([0, 6])
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Error (mV)')

    ax1 = plt.subplot(222)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    startT = 0
    endT = 2000
    ax1.plot(plotT[startT:endT], neuronVArray[startT:endT], alpha=.7, label='$NEURON$', color=viridis_cmap(1))
    ax1.plot(plotT[startT:endT], spinnaker[startT:endT], alpha=.7, label='$SpiNNaker$', color=viridis_cmap(0.6))
    ax1.plot(plotT[startT:endT], jib2_array[startT:endT], alpha=.7, label='$Jib2$', color=viridis_cmap(0.2))
    ax1.set_yticks([-75, 25])
    ax1.set_xticks([0, 200])

    ax1 = plt.subplot(224)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.plot(plotT[startT:endT], error_jib2[startT:endT], color=viridis_cmap(0.2), alpha=0.7, label = '$Jib2$')
    ax1.plot(plotT[startT:endT], error_spiNNaker[startT:endT], color=viridis_cmap(0.6), alpha=0.7, label = '$SpiNNaker$')
    ax1.set_xticks([0, 200])
    ax1.set_xlabel('Time (ms)')
    plt.show()
    return


def main():
    comp_neuron_with_spin()


if __name__ == '__main__':
    main()
