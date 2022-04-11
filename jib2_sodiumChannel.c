#include <qpe.h>  // required qpe support functions
#include <nmu.h>  // neuromorphic accelerators


int main ()
{
  nmu_exp_log_control(14);
 // variables to be uploaded MUST be static
  //NOTE: otherwise their value can change before retrieval
  static char *   msg = "float exp(x) results:";
  static uint32_t num_data;
  static float    data[1000];

  // initialise variables and constants
  # define constants
  static float cm = 1;
  static float dt = 0.1;
  static float gNa = 0.12;
  static float gL = 0.0003;
  static float El = -54.3;
  static float ENa = 50;
  static float m = 0;
  static float h = 0;
  static float i_offset = 5;
  static float V_soma = -65;
  static float V_inf;
  static float tau_V;
  static float alpha_n;
  static float beta_n;
  static float tau_m;
  static float m_inf;
  static float tau_h;
  static float h_inf;
  static float totGi;
  static float totGiE;
  uint32_t count_refrac = 0;
  static float v_thresh = -50;
  static float v_reset = -65;
  static float tau_refrac = 1;

  //alpha_n = (0.01*(V_soma+55)) / (1 - hwexpf(-0.1*(V_soma+55)));
  //beta_n = 0.125*hwexpf(-0.0125*(V_soma+65));
  //tau_n = 1/(alpha_n+beta_n);
  //n_inf = alpha_n/(alpha_n+beta_n);
  n = 0;


  num_data = 200;
  //run through timesteps`
  for (uint32_t i = 0; i < num_data; i++) {
      // record value of V
      data[i] = gNa*m*m*m*h*(V_soma-ENa);
    
    	alpha_m = (0.1*(V_soma+40)) / (1 - hwexpf(-0.1*(V_soma+40)))
      beta_m = 4*hwexpf(-0.0556*(V_soma+65))
      tau_m = 1/(alpha_m+beta_m)
      m_inf = alpha_m/(alpha_m+beta_m)
      m = m_inf + (m-m_inf)*hwexpf(-dt/tau_m)

      alpha_h = 0.07*hwexpf(-0.05*(V_soma+65))
      beta_h = 1 / (1+hwexpf(-0.1*(V_soma+35)))
      tau_h = 1/(alpha_h+beta_h)
      h_inf = alpha_h/(alpha_h+beta_h)
      h = h_inf + (h-h_inf)*hwexpf(-dt/tau_h)

      if (count_refrac == 0){
      totGi = gL + gNa*(m*m*m*h);
      totGiE = gL*El + gNa*(m*m*m*h)*ENa;

      V_inf = (totGiE + i_offset) / totGi;
      tau_V = cm/totGi;
      V_soma = V_inf + (V_soma - V_inf)*hwexpf(-dt/tau_V);
      }

      else {
      count_refrac -= 1;
      }

      if (V_soma > v_thresh){
      V_soma = v_reset;
      count_refrac = tau_refrac;
      }


  }

// pass variable addresses as results
//NOTE: arrays and strings are already addresses!
exit_res[0] = (void *) msg;
exit_res[1] = (void *) &num_data;  // scalars need the & operator
exit_res[2] = (void *) data;

// send exit code to host - can be ignored
//NOTE: by convention, 0 indicates successfull finish
return 0;

// compute exp(x) using the hardware accelerator
//  num_data = 4;
//  for (uint32_t i = 0; i < num_data; i++) {
//    data[i] = hwexpf (1.0 * i);
//  }
