  // define constants and variables
  s1615 dt = 0.1;
  s1615 V_dend;
  V_dend = -65;
  s1615 Ie_dend = 35.80986219567645;
  s1615 dCaAP_count = 0;
  s1615 t_dCaAP;
  t_dCaAP = -200;
  s1615 K = 0;
  s1615 i_dCaAP = 0;
  s1615 v_th = -36;
  s1615 refract_period = 200;
  s1615 AB_refract_period = 77.7;
  s1615 denom = -0.08547009;
  s1615 w = 0.2;
  s1615 this_t = -0.05;
  int count = 0;
  s1615 AB = 0;
  s1615 totGE_dend = -70;
  s1615 C_dend = -1;
  s1615 F_dend;
  s1615 b_dend = 0;
  s1615 c_dend = -0.1;
  s1615 f_dend;
  s1615 c_x_dend;
  s1615 f_x_dend;
  s1615 dV_dend;
  s1615 this_Ie_dend;

  s1615 V_soma;
  s1615 v_rest = -100;
  V_soma = -65;
  s1615 m = 0;
  s1615 n = 0;
  s1615 h = 0;
  s1615 totG_soma;
  s1615 totGE_soma;
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
  s1615 g_soma = 0.1;
  s1615 C_soma;
  s1615 D_soma = 0.1;
  s1615 F_soma;
  s1615 c_soma;
  s1615 d_soma = 0.01;
  s1615 f_soma;
  s1615 c_x_soma;
  s1615 f_x_soma;
  s1615 dV_soma;


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
  
  
  num_data = 10000;
  
  for (int i = 0; i < num_data; i++) {

      // use values of m, n and h to calculate current and update membrane voltage
      totG_soma =  (gL + gK * (n*n*n*n) + gNa * (m*m*m) * h)*1000;
      totGE_soma = (gL * El + gK * (n*n*n*n) * Ek + gNa * (m*m*m) * h * ENa)*1000;

      // check whether a dCaAP has been fired
      this_t = (i*dt) - (dt/2);


		  A_dCaAP = A_dCaAP + (1. - expk(dt * (((((1.0 - A_dCaAP)) + (A_dCaAP) * (((- 1.0))))) / tau_A))) * (- A_dCaAP);
// 	   B_dCaAP = B_dCaAP * (expk(x2 - y2*B_dCaAP));

		  B_dCaAP = B_dCaAP + (1. - expk(dt * ((((1.0) * ((1.0 - B_dCaAP)) + (B_dCaAP) * (((- 1.0))))) / tau_B))) * (- B_dCaAP);


		  if (B_dCaAP == 0 && this_t > t_dCaAP + sigma_diff)	{
			  B_dCaAP = 0.001;
		  }

		  if (this_t > t_dCaAP + refract_period && V_dend > v_th) {
			  t_dCaAP = this_t;
			  A_dCaAP = 0.001k;
			  B_dCaAP = 0;
        K = expk((V_dend - v_th) * denom);

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
}
