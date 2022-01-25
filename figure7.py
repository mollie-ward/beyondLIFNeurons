import numpy as np
import matplotlib.pyplot as plt
import functions as functions
import matplotlib as mlib

# fig1 = plt.figure(facecolor='white')
viridis_cmap = mlib.cm.get_cmap('viridis')
purd_cmap = mlib.cm.get_cmap('PuRd')
mlib.rcParams.update({'font.size': 24})


def sim_soma_dcaap():
    these_synapses = [10, 20, 40, 80]
    count = 1

    for synapses in these_synapses:
        num_synapses = synapses
        total_time = 2
        firing_frequency = 20
        spike_time_array = np.zeros([num_synapses, firing_frequency])

        for synapse in range(num_synapses):
            spike_time_array[synapse, :] = np.sort(np.random.randint(0, total_time * 10000 + 1, firing_frequency))

        synapse_ps = np.zeros(num_synapses)

        num_compartments = 2
        v = -65
        z = 1

        g = np.empty((num_compartments, 2))
        g[0, 0] = 0
        g[0, 1] = 0.1
        g[1, 0] = 0
        g[1, 1] = 0


        dt = 0.1


        time = total_time * 1000
        # ms of simulation
        t = time / dt
        Ie = np.zeros((num_compartments, int(t)+1))

        soma_spike_count = 0

        Ie[1, :] = 0

        # start dendrite off
        V_dend = -65
        dCaAP_count = 0
        A_dCaAP = 0
        B_dCaAP = 0
        t_dCaAP = -200
        K = 0
        A_dCaAP, B_dCaAP, i, dCaAP_count, t_dCaAP, K = dCaAP(V_dend, t_dCaAP, dCaAP_count, A_dCaAP, B_dCaAP, -dt/2, dt, K)

        # start soma off
        V_soma = -65
        m = 0.0529
        n = 0.3177
        h = 0.5961
        totG_soma, totGE_soma, m, n, h = gating_variables_hh(V_soma, m, n, h, dt)

        synaptic_current = []
        ca_current = []
        voltage_array = []
        dend_voltage_array = []
        voltage_array.append(V_soma)
        dend_voltage_array.append(V_dend)
        timesteps = np.linspace(0, total_time*1000, total_time*10000)

        for timestep in range(1, int(t)):
            t = round(timesteps[timestep], 5)
            old_V_soma = V_soma
            voltage_array.append(V_soma)
            dend_voltage_array.append(V_dend)

            # update synaptic current
            synapse_ps, i_synapses = functions.update_synapse(synapse_ps, timestep, spike_time_array, V_dend)
            synaptic_current.append(i_synapses)
            Ie[1, timestep] += i_synapses*100
            
            # update dendritic current
            totG_dend, totGE_dend, i = gating_variables_passive(V_dend)
            A_dCaAP, B_dCaAP, i_dCaAP, dCaAP_count, t_dCaAP, K = dCaAP(V_dend, t_dCaAP, dCaAP_count, A_dCaAP, B_dCaAP, t-(dt/2), dt, K)
            ca_current.append(i_dCaAP)
            Ie[1, timestep] -= i_dCaAP*100

            # update somatic current
            totG_soma, totGE_soma, m, n, h = gating_variables_hh(V_soma, m, n, h, dt)

            # update voltages
            V_dend, V_soma = voltage_update_twoComp(V_dend, V_soma, dt, totG_dend, totGE_dend, totG_soma, totGE_soma, g, Ie[:, timestep])

            if V_soma > 20 and old_V_soma < 20:
                soma_spike_count += 1


        ax1 = plt.subplot(2, 4, count)
        ax1.plot(timesteps[1000:], dend_voltage_array[1000:], color=purd_cmap(0.85), alpha=1, label='$V_{dend}$', linewidth=1.5)
        ax1.axes.get_yaxis().set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.yaxis.set_ticks_position('left')
        ax1.set_ylim([-80, 50])
        ax1.set_xticks([0, 2000])
        ax1.set_xticklabels(['0', '2'])
        
        if count == 1:
            ax1.axes.get_yaxis().set_visible(True)
            ax1.set_yticks([-50, 50])
            ax1.set_xlabel("Time (s)")
            ax1.set_ylabel("Membrane voltage (mV)")

        if count == 4:
            plt.legend(frameon=False)

        ax1 = plt.subplot(2, 4, count+4)
        ax1.plot(timesteps[1000:], voltage_array[1000:], color=viridis_cmap(0.15), alpha=1, label='$V_{soma}$', linewidth=1.5)
        ax1.axes.get_yaxis().set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.set_ylim([-80, 50])
        ax1.set_xticks([0, 2000])
        ax1.set_xticklabels(['0', '2'])

        if count == 1:
            ax1.axes.get_yaxis().set_visible(True)       
            ax1.set_yticks([-50, 50])
            ax1.set_xlabel("Time (s)")
            ax1.set_ylabel("Membrane voltage (mV)")
        if count == 4:
            plt.legend(frameon=False)

        count += 1


    plt.show()

    return
  

def dCaAP(v, t_dCaAP, dCaAP_count, A, B, t, dt, K):
    v_th = -36
    refract_period = 200
    tauA = 3
    tauB = 0.4
    D = 0.3
    sigma_diff = 21
    v_rest = -75
    denom = -1 / ((v_th - v_rest)*D)
    w = 3# change to 3
    x = dt * (1/tauA)
    y = 2 * x

    x2 = dt * (1/tauB)
    y2 = 2 * x2

    A = A * np.exp(x - y*A)
    B = B * (np.exp(x2 - y2*B))

    if B == 0 and t > t_dCaAP + sigma_diff and w > 0:
        B = 0.001

    if t > t_dCaAP + refract_period and v > v_th and w > 0:
        t_dCaAP = t
        A = 0.001
        B = 0
        K = np.exp((v - v_th) * denom)
        if K > 1:
            K = 1
        dCaAP_count += 1
        
    i = -(A - B) * w * K

    return A, B, i, dCaAP_count, t_dCaAP, K
 


def gating_variables_hh(V, m, n, h, dt):

    gNa = 0.12
    gK = 0.036
    gL = 0.0003
    El = -54.3
    ENa = 50
    Ek = -77

    # for compartment in range(0, len(V)):
    alpha_m = (0.1 * (V + 40)) / (1 - np.exp(-0.1 * (V + 40)))
    beta_m = 4 * np.exp(-0.0556 * (V + 65))
    tau_m = 1 / (alpha_m + beta_m)
    m_inf = alpha_m / (alpha_m + beta_m)
    m = m_inf + (m - m_inf) * np.exp(-dt / tau_m)

    alpha_h = 0.07 * np.exp(-0.05 * (V + 65))
    beta_h = 1 / (1 + np.exp(-0.1 * (V + 35)))
    tau_h = 1 / (alpha_h + beta_h)
    h_inf = alpha_h / (alpha_h + beta_h)
    h = h_inf + (h - h_inf) * np.exp(-dt / tau_h)

    alpha_n = (0.01 * (V + 55)) / (1 - np.exp(-0.1 * (V + 55)))
    beta_n = 0.125 * np.exp(-0.0125 * (V + 65))
    tau_n = 1 / (alpha_n + beta_n)
    n_inf = alpha_n / (alpha_n + beta_n)
    n = n_inf + (n - n_inf) * np.exp(-dt / tau_n)

    totGi = gL + gK*(n**4) + gNa*(m**3)*h
    totGiE = gL*El + gK*(n**4)*Ek + gNa*(m**3)*h*ENa

    return totGi*1000, totGiE*1000, m, n, h
  
  
def gating_variables_passive(V):
    gL_bar = 0.001
    eL = -70
    i = gL_bar*(V - eL)

    totG = gL_bar
    totGE = gL_bar*eL

    return totG*1000, totGE*1000, i


def voltage_update_twoComp(V_dend, V_soma, dt, totG_dend, totGE_dend, totG_soma, totGE_soma, g, Ie):
    cm = 1
    # update somatic
    B_soma = 0
    C_soma = -(1 / cm) * (totG_soma + g[0, 1])
    D_soma = (1/cm) * g[0, 1]
    F_soma = (1/cm) * (totGE_soma + Ie[0])

    b_soma = 0
    c_soma = C_soma * dt
    d_soma = D_soma * dt
    f_soma = (F_soma + C_soma*V_soma + D_soma*V_dend) * dt

    # update dendrite
    B_dend = (1/cm) * g[1, 0]
    C_dend = -(1 / cm) * (totG_dend + g[1, 0])
    D_dend = 0
    F_dend = (1/cm) * (totGE_dend + Ie[1])

    b_dend = B_dend * dt
    c_dend = C_dend * dt
    d_dend = 0
    f_dend = (F_dend + B_dend*V_soma + C_dend*V_dend) * dt

    c_x_soma = c_soma
    f_x_soma = f_soma
    c_x_dend = c_dend + (b_dend*d_soma)/(1-c_x_soma)
    f_x_dend = f_dend + (b_dend*f_soma)/(1-c_x_soma)

    dV_dend = f_x_dend / (1 - c_x_dend)
    dV_soma = (d_soma*dV_dend + f_x_soma) / (1 - c_x_soma)

    V_dend += dV_dend
    V_soma += dV_soma

    return V_dend, V_soma
  
  
def update_synapse(ps_vals, timestep, spike_times, V):

    #update synaptic current
    Gs = 0.002

    # num_synapses = 35
    exp_decay = np.exp(-0.1/12)
    # test_exp = np.exp(-0.1/1.8) - np.exp(-0.1/0.3)
    current = np.empty(np.size(ps_vals))

    for synapse in range(np.size(ps_vals)):
        ps_vals[synapse] *= exp_decay
        these_spike_times = spike_times[synapse,:]

        if timestep in these_spike_times:
            ps_vals[synapse] += (1 - ps_vals[synapse])

        current[synapse] = Gs * ps_vals[synapse] * (0 - V)

    total_current = np.sum(current)

    return ps_vals, total_current
  

def main():
    sim_soma_dcaap()


if __name__ == '__main__':
    main()
