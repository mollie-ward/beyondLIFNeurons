#include <qpe.h>  // required qpe support functions
#include <nmu.h>  // neuromorphic accelerators

int main ()
{
  static uint32_t num_data;
  static uint32_t num_data;

  // define constants and variables
  static float dt = 0.1;
  static float V_dend;
  V_dend = -65;
  static float Ie_dend = 35.80986219567645;
  static float dCaAP_count = 0;
  static float t_dCaAP;
  t_dCaAP = -200;
  static float K = 0;
  static float i_dCaAP = 0;
  static float v_th = -36;
  static float refract_period = 200;
  static float AB_refract_period = 77.7;
  static float denom = -0.08547009;
  static float w = 0.2;
  static float this_t = -0.05;
  static uint32_t count = 0;
  static float AB = 0;
  static float totGE_dend = -70;
  static float C_dend = -1;
  static float F_dend;
  static float b_dend = 0;
  static float c_dend = -0.1;
  static float f_dend;
  static float c_x_dend;
  static float f_x_dend;
  static float dV_dend;
  static float this_Ie_dend;

  static float V_soma;
  static float v_rest = -100;
  V_soma = -65;
  static float m = 0;
  static float n = 0;
  static float h = 0;
  static float totG_soma;
  static float totGE_soma;
  static float gNa = 0.12;
  static float gK = 0.036;
  static float gL = 0.0003;
  static float El = -54.3;
  static float ENa = 50;
  static float Ek = -77;
  static float alpha_m;
  static float beta_m;
  static float alpha_h;
  static float beta_h;
  static float alpha_n;
  static float beta_n;
  static float tau_m;
  static float tau_h;
  static float tau_n;
  static float m_inf;
  static float h_inf;
  static float n_inf;
  static float Ie_soma = 7.957747154594768;
  static float g_soma = 0.1;
  static float C_soma;
  static float D_soma = 0.1;
  static float F_soma;
  static float c_soma;
  static float d_soma = 0.01;
  static float f_soma;
  static float c_x_soma;
  static float f_x_soma;
  static float dV_soma;


  
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

  alpha_n = (0.01*(V_soma+55)) / (1 - hwexpf(-0.1*(V_soma+55)));
  beta_n = 0.125*hwexpf(-0.0125*(V_soma+65));
  tau_n = 1/(alpha_n+beta_n);
  n_inf = alpha_n/(alpha_n+beta_n);
  n = n_inf;
  
  
  num_data = 10000;
  
  for (uint32_t i = 0; i < num_data; i++) {

      // use values of m, n and h to calculate current and update membrane voltage
      totG_soma =  (gL + gK * (n*n*n*n) + gNa * (m*m*m) * h)*1000;
      totGE_soma = (gL * El + gK * (n*n*n*n) * Ek + gNa * (m*m*m) * h * ENa)*1000;

      // check whether a dCaAP has been fired
      this_t = (i*dt) - (dt/2);


		  A_dCaAP = A_dCaAP + (1. - hwexpf(dt * (((((1.0 - A_dCaAP)) + (A_dCaAP) * (((- 1.0))))) / tau_A))) * (- A_dCaAP);
// 	   B_dCaAP = B_dCaAP * (hwexpf(x2 - y2*B_dCaAP));

		  B_dCaAP = B_dCaAP + (1. - hwexpf(dt * ((((1.0) * ((1.0 - B_dCaAP)) + (B_dCaAP) * (((- 1.0))))) / tau_B))) * (- B_dCaAP);


		  if (B_dCaAP == 0 && this_t > t_dCaAP + sigma_diff)	{
			  B_dCaAP = 0.001;
		  }

		  if (this_t > t_dCaAP + refract_period && V_dend > v_th) {
			  t_dCaAP = this_t;
			  A_dCaAP = 0.001k;
			  B_dCaAP = 0;
        K = hwexpf((V_dend - v_th) * denom);

        if (K > 1) {
          K = 1;
        }

        dCaAP_count += 1;
		}

		  i_dCaAP = -(A_dCaAP - B_dCaAP) * w * K;
      this_Ie_dend = Ie_dend - i_dCaAP*100;

    
      C_soma = -(totG_soma + g_soma);
      F_soma = totGE_soma + Ie_soma;

      c_soma = C_soma*dt;
      f_soma = (F_soma + C_soma*V_soma + D_soma*V_dend) * dt;

      F_dend = totGE_dend + this_Ie_dend;
      f_dend = (F_dend + C_dend*V_dend) * dt;

      c_x_soma = c_soma;
      f_x_soma = f_soma;

      c_x_dend = c_dend + (b_dend*d_soma)/(1 - c_x_soma);
      f_x_dend = f_dend + (b_dend*f_soma)/(1 - c_x_soma);

      dV_dend = f_x_dend / (1 - c_x_dend);
      dV_soma = (d_soma*dV_dend + f_x_soma) / (1 - c_x_soma);

      V_dend += dV_dend;
      V_soma += dV_soma;

      // update m, n and h using updated membrane voltage 
      // update m, n and h using updated membrane voltage 
      alpha_m = (0.1*(V_soma+40)) / (1 - hwexpf(-0.1*(V_soma+40)));
      beta_m =  4*hwexpf(-0.0556*(V_soma+65));
      tau_m = 1/(alpha_m+beta_m);
      m_inf = alpha_m/(alpha_m+beta_m);
      m = m_inf + (m-m_inf)*hwexpf(-dt/tau_m);


      alpha_h = 0.07*hwexpf(-0.05*(V_soma+65));
      beta_h =  1 / (1+hwexpf(-0.1*(V_soma+35)));
      tau_h = 1/(alpha_h+beta_h);
      h_inf = alpha_h/(alpha_h+beta_h);
      h = h_inf + (h-h_inf)*hwexpf(-dt/tau_h);

      alpha_n = (0.01*(V_soma+55)) / (1 - hwexpf(-0.1*(V_soma+55)));
      beta_n = 0.125*hwexpf(-0.0125*(V_soma+65));
      tau_n = 1/(alpha_n+beta_n);
      n_inf = alpha_n/(alpha_n+beta_n);
      n = n_inf + (n-n_inf)*hwexpf(-dt/tau_n);
    
  }
}
