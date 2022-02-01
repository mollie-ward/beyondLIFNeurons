// initialise variables and constants
s1615 dt = 0.1;
s1615 V_soma;
V_soma = -65;
s1615 m = 0;
s1615 n = 0;
s1615 h = 0;
s1615 totG_soma = 0;
s1615 totGE_soma = 0;
s1615 gNa = 0.12;
s1615 gK = 0.036;
s1615 gL = 0.0003;
s1615 El = -54.3;
s1615 ENa = 50;
s1615 Ek = -77;
s1615 alpha_m;
s1615 beta_m;
s1615 alpha_h;
s1615 beta_h;
s1615 alpha_n;
s1615 beta_n;
s1615 tau_m;
s1615 tau_h;
s1615 tau_n;
s1615 m_inf;
s1615 h_inf;
s1615 n_inf;
s1615 Ie_soma = 7.957747154594768;
s1615 C_soma;
s1615 F_soma;
s1615 c_soma;
s1615 f_soma;
s1615 dV_soma;

// initialse values of m, n and h as m_inf, n_inf and h_inf
alpha_m = (0.1*(V_soma+40)) / (1 - expk(-0.1*(V_soma+40)));
beta_m =  4*expk(-0.0556*(V_soma+65));
tau_m = 1/(alpha_m+beta_m);
m_inf = alpha_m/(alpha_m+beta_m);
m = m_inf;

alpha_h = 0.07*expk(-0.05*(V_soma+65));
beta_h =  1 / (1+expk(-0.1*(V_soma+35)));
tau_h = 1/(alpha_h+beta_h);
h_inf = alpha_h/(alpha_h+beta_h);
h = h_inf;

alpha_n = (0.01*(V_soma+55)) / (1 - expk(-0.1*(V_soma+55)));
beta_n = 0.125*expk(-0.0125*(V_soma+65));
tau_n = 1/(alpha_n+beta_n);
n_inf = alpha_n/(alpha_n+beta_n);
n = n_inf;


num_data = 20001;
//run through timesteps`
for (int i = 0; i < num_data; i++) {
    // record value of V
    data[i] = V_soma;

    // use values of m, n and h to calculate current and update membrane voltage
    totG_soma =  (gL + gK * (n*n*n*n) + gNa * (m*m*m) * h)*1000;
    totGE_soma = (gL * El + gK * (n*n*n*n) * Ek + gNa * (m*m*m) * h * ENa)*1000;
    C_soma = -1*totG_soma;
    F_soma = totGE_soma + Ie_soma;
    c_soma = C_soma*dt;
    f_soma = (F_soma + C_soma*V_soma)*dt;
    c_soma = 1-c_soma;
    dV_soma = f_soma/c_soma;
    V_soma += dV_soma;

    // update m, n and h using updated membrane voltage 
    alpha_m = (0.1*(V_soma+40)) / (1 - expk(-0.1*(V_soma+40)));
    beta_m =  4*expk(-0.0556*(V_soma+65));
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    m = m_inf + (m-m_inf)*expk(-dt/tau_m);


    alpha_h = 0.07*expk(-0.05*(V_soma+65));
    beta_h =  1 / (1+expk(-0.1*(V_soma+35)));
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    h = h_inf + (h-h_inf)*expk(-dt/tau_h);

    alpha_n = (0.01*(V_soma+55)) / (1 - expk(-0.1*(V_soma+55)));
    beta_n = 0.125*expk(-0.0125*(V_soma+65));
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    n = n_inf + (n-n_inf)*expk(-dt/tau_n);
    }
