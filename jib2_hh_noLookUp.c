#include <qpe.h>  // required qpe support functions
#include <nmu.h>  // neuromorphic accelerators

int main ()
{     
  // variables to be uploaded MUST be static
  //NOTE: otherwise their value can change before retrieval
  static char *   msg = "float exp(x) results:";
  static uint32_t num_data;
  static float    data[20001];

  // initialise variables and constants
  static float dt = 0.1;
  static float V_soma;
  V_soma = -65;
  static float m = 0;
  static float n = 0;
  static float h = 0;
  static float totG_soma = 0;
  static float totGE_soma = 0;
  static float gNa = 0.12;
  static float gK = 0.036;
  static float gL = 0.0003;
  static float El = -54.3;
  static float ENa = 50;
  static float Ek = -77;
  static float tau_m;
  static float tau_h;
  static float tau_n;
  static float m_inf;
  static float h_inf;
  static float n_inf;
  static float Ie_soma = 7.957747154594768;
  static float C_soma;
  static float F_soma;
  static float c_soma;
  static float f_soma;
  static float dV_soma;

  // initialse values of m, n and h as m_inf, n_inf and h_inf
  alpha_m = (0.1*(V_soma+40)) / (1 - hwexpf(-0.1*(V_soma+40)));
  beta_m =  4*hwexpf(-0.0556*(V_soma+65));
  tau_m = 1/(alpha_m+beta_m);
  m_inf = alpha_m/(alpha_m+beta_m);
  m = m_inf;

  alpha_h = 0.07*hwexpf(-0.05*(V_soma+65));
  beta_h =  1 / (1+hwexpf(-0.1*(V_soma+35)));
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
  for (uint32_t i = 0; i < num_data; i++) {
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
    
      n_pow = gK*n*n*n*n;
      m_h_pow = gNa*m*m*m*h;
      totG_soma = (gL + n_pow + m_h_pow)*1000;
      totGE_soma = (gL*El + n_pow*Ek + m_h_pow*ENa)*1000;
      }

// pass variable addresses as results
//NOTE: arrays and strings are already addresses!
exit_res[0] = (void *) msg;
exit_res[1] = (void *) &num_data;  // scalars need the & operator
exit_res[2] = (void *) data;

// send exit code to host - can be ignored
//NOTE: by convention, 0 indicates successfull finish
return 0;
}
