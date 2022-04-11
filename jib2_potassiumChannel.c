#include <qpe.h>  // required qpe support functions
#include <nmu.h>  // neuromorphic accelerators


int main ()
{
  nmu_exp_log_control(14);
 // variables to be uploaded MUST be static
  //NOTE: otherwise their value can change before retrieval
  static char *   msg = "float exp(x) results:";
  static uint32_t num_data;
  static float    data[200];

  // initialise variables and constants
  # define constants
  static float cm = 1;
  static float dt = 0.1;
  static float gK = 0.036;
  static float gL = 0.0003;
  static float El = -54.3;
  static float Ek = -77;
  static float n = 0.32;
  static float i_offset = 5;
  static float V_soma = -65;
  static float V_inf;
  static float tau_V;
  static float alpha_n;
  static float beta_n;
  static float tau_n;
  static float n_inf;
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
      data[i] = gK*n*n*n*n*(V_soma-Ek);

      alpha_n = (0.01*(V_soma+55)) / (1 - hwexpf(-0.1*(V_soma+55)));
      beta_n = 0.125*hwexpf(-0.0125*(V_soma+65));
      tau_n = 1/(alpha_n+beta_n);
      n_inf = alpha_n/(alpha_n+beta_n);
      n = n_inf + (n-n_inf)*hwexpf(-dt/tau_n);

      if (count_refrac == 0){
      totGi = gL + gK*(n*n*n*n);
      totGiE = gL*El + gK*(n*n*n*n)*Ek;

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
