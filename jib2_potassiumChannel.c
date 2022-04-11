#include <qpe.h>  // required qpe support functions
#include <nmu.h>  // neuromorphic accelerators

int main ()
{     
  // variables to be uploaded MUST be static
  //NOTE: otherwise their value can change before retrieval
  static char *   msg = "float exp(x) results:";
  static uint32_t num_data;
  static float    data[10001];

  // initialise variables and constants
  # define constants
  static float cm = 1;
  static float dt = 0.1;
  static float gNa = 0.12;
  static float gK = 0.036;
  static float gL = 0.0003;
  static float El = -54.3;
  static float ENa = 50;
  static float Ek = -77;
  static float n = 0.32;
  static float i_offset = 10;
  static float V_soma = -65;
  static float V_inf;
  static float tau_V;


  alpha_n = (0.01*(V_soma+55)) / (1 - hwexpf(-0.1*(V_soma+55)));
  beta_n = 0.125*hwexpf(-0.0125*(V_soma+65));
  tau_n = 1/(alpha_n+beta_n);
  n_inf = alpha_n/(alpha_n+beta_n);
  n = n_inf;
  

  num_data = 10001;
  //run through timesteps`
  for (uint32_t i = 0; i < num_data; i++) {
      // record value of V
      data[i] = V_soma;

      alpha_n = (0.01*(V_soma+55)) / (1 - hwexpf(-0.1*(V_soma+55)));
      beta_n = 0.125*hwexpf(-0.0125*(V_soma+65));
      tau_n = 1/(alpha_n+beta_n);
      n_inf = alpha_n/(alpha_n+beta_n);
      n = n_inf + (n-n_inf)*hwexpf(-dt/tau_n);
    
      totGi = gL + gK*(n**4);
      totGiE = gL*El + gK*(n**4)*Ek;

      V_inf = (totGiE + i_offset) / totGi;
      tau_V = cm/totGi;    
    
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


