#!/usr/bin/python
# Here we make a direct comparison between NEURON and python for a simple one compartment model

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlib

# fig1 = plt.figure(facecolor='white')
viridis_cmap = mlib.cm.get_cmap('viridis')
# viridis_cmap = mlib.cm.get_cmap('plasma')
purd_cmap = mlib.cm.get_cmap('PuRd')


mlib.rcParams.update({'font.size': 24})

def python_sim_hh(inject, length, time):
    voltage = []
    m_array = []
    n_array = []
    h_array = []

    soma_v = -65
    dt = 0.1
    time = time  # ms of simulation
    t = time / dt
    ie = np.zeros(int(t) + 1)
    m = 0
    h = 0
    n = 0
    table = lookup_table_hh()

    # initialise m, n and h as minf, ninf, hinf from LUT
    # calculate remainder for interpolation
    rem = (soma_v - -100) - int(soma_v - -100)
    # calculate indies for LUT
    indices = [int(soma_v - -100), int(soma_v - -100) + 1]
    # calculate minf, hinf, ninf
    m_inf = np.diff(table[2, indices]) * rem + table[2, indices[0]]
    h_inf = np.diff(table[4, indices]) * rem + table[4, indices[0]]
    n_inf = np.diff(table[6, indices]) * rem + table[6, indices[0]]
    # in INITIAL block in NEURON code m=minf, h=hinf, n=ninf
    m = m_inf[0]
    h = h_inf[0]
    n = n_inf[0]


    # convert the Ie in NEURON to Ie/A - nA to mA/cm2
    # Ie = 0.9
    # # ie[4000:8000] = (Ie/A)
    # area = 3.14159265358978952e-6
    # Ie_mA = 9e-7/area
    Ie = inject/1256.6370614359173*100000
    ie[0:length+1] = Ie

    for timestep in range(0, int(t) + 1):
        # store values for plotting
        voltage.append(soma_v)
        m_array.append(m)
        n_array.append(n)
        h_array.append(h)

        # update ion channel currents
        totG, totGE = BREAKPOINT(m, h, n, soma_v)

        # update voltage

        C = - totG
        F = totGE + ie[timestep]
        c = C * dt
        f = (F + C * soma_v) * dt
        c_x = c
        f_x = f
        dV = f_x / (1 - c_x)
        soma_v = soma_v + dV

        # update ion channel params based on new voltage
        totG, totGE, m, n, h = hh_lookup(soma_v, m, n, h, table)

    return voltage, m_array, n_array, h_array


def soma_figure():
    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    time = 40
    single_spike, m, n, h = python_sim_hh(0.03769911184, time*10, time)
    timesteps = np.linspace(0, int(time+1), int(time * 10)+1)
    ax1.plot(timesteps[:], single_spike[:], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth = 1.5)
    ax1.set_yticks([-75, 20])
    ax1.set_ylabel('Membrane voltage (mV)')

    ax1 = plt.subplot(223)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    time = 40

    timesteps = np.linspace(0, int(time+1), int(time * 10)+1)
    ax1.plot(timesteps[:], m[:], alpha=1, label='$m$', linewidth = 1.5)
    ax1.plot(timesteps[:], n[:], alpha=1, label='$n$', linewidth = 1.5)
    ax1.plot(timesteps[:], h[:], alpha=1, label='$h$', linewidth = 1.5)
    ax1.set_yticks([0, 1])
    ax1.legend(frameon=False)
    ax1.set_ylabel('Gating variables')
    ax1.set_xlabel('Time (ms)')

    ax3 = plt.subplot(222)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.yaxis.set_ticks_position('left')
    small_input, m, n, h = python_sim_hh(0.3769911184, time*10, time)
    ax3.plot(timesteps[:], small_input[:], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth = 1.5)
    ax3.set_yticks([-75, 20])
    ax3.set_xlabel('Time (ms)')

    ax1 = plt.subplot(224)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    time = 40

    timesteps = np.linspace(0, int(time+1), int(time * 10)+1)
    ax1.plot(timesteps[:], m[:], alpha=1, label='$m$', linewidth = 1.5)
    ax1.plot(timesteps[:], n[:], alpha=1, label='$n$', linewidth = 1.5)
    ax1.plot(timesteps[:], h[:], alpha=1, label='$h$', linewidth = 1.5)
    ax1.set_yticks([0, 1])
    ax1.set_xlabel('Time (ms)')


def lookup_table_hh():
    v_values = np.arange(-100, 100, 1)
    table = np.empty([7, len(v_values)])
    table[0, :] = v_values
    count = 0

    for vt in v_values:
        alpha_m = .1 * vtrap(-(vt + 40), 10)
        beta_m = 4 * np.exp(-(vt + 65) / 18)
        tau_m = 1 / (alpha_m + beta_m)
        m_inf = alpha_m / (alpha_m + beta_m)

        alpha_h = 0.07 * np.exp(-(vt+65)/20)
        beta_h = 1 / (np.exp(-(vt+35)/10) + 1)
        tau_h = 1 / (alpha_h + beta_h)
        h_inf = alpha_h / (alpha_h + beta_h)

        alpha_n = .01*vtrap(-(vt+55), 10)
        beta_n = .125*np.exp(-(vt+65)/80)
        tau_n = 1 / (alpha_n + beta_n)
        n_inf = alpha_n / (alpha_n + beta_n)

        table[1, count] = tau_m
        table[2, count] = m_inf
        table[3, count] = tau_h
        table[4, count] = h_inf
        table[5, count] = tau_n
        table[6, count] = n_inf

        count = count + 1

    return table


def hh_lookup(v, m, n, h, table):
    gnabar = 0.12
    gkbar = 0.036
    gl = 0.0003
    el = -54.3
    ena = 50
    ek = -77


    rem = (v - -100) - int(v - -100)
    indices = [int(v - -100), int(v - -100)+1]

    tau_m = np.diff(table[1, indices]) * rem + table[1, indices[0]]
    m_inf = np.diff(table[2, indices]) * rem + table[2, indices[0]]
    tau_h = np.diff(table[3, indices]) * rem + table[3, indices[0]]
    h_inf = np.diff(table[4, indices]) * rem + table[4, indices[0]]
    tau_n = np.diff(table[5, indices]) * rem + table[5, indices[0]]
    n_inf = np.diff(table[6, indices]) * rem + table[6, indices[0]]

    m = m_inf[0] - (-m + m_inf[0]) * np.exp(-0.1 / tau_m[0])
    n = n_inf[0] - (-n + n_inf[0]) * np.exp(-0.1 / tau_n[0])
    h = h_inf[0] - (-h + h_inf[0]) * np.exp(-0.1 / tau_h[0])

    totG = gl + gkbar * (n ** 4) + gnabar * (m ** 3) * h
    totGE = gl * el + gkbar * (n ** 4) * ek + gnabar * (m ** 3) * h * ena

    return totG*1000, totGE*1000, m, n, h



def BREAKPOINT(m,h,n,v):
    gnabar = 0.12
    gkbar = 0.036
    gl = 0.0003
    el = -54.3
    ena = 50
    ek = -77

    totG = gl + gkbar * (n ** 4) + gnabar * (m ** 3) * h
    totGE = gl * el + gkbar * (n ** 4) * ek + gnabar * (m ** 3) * h * ena


    return totG*1000, totGE*1000


def vtrap(x,y):
    # Traps for 0 in denominator of rate eqns.
    if abs(x / y) < 1e-6:
        vtrap = y*(1 - x/y/2)
    else:
        vtrap = x/(np.exp(x/y) - 1)
    return vtrap


def main():
    soma_figure()
    plt.show()


if __name__ == '__main__':
    main()
