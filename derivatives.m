% derivatives.m
% This function calculates the derivatives for all state variables in the rod photoreceptor model.
% It is designed to be called by a MATLAB ODE solver (e.g., ode15s).

function dYdt = derivatives(t, Y, params)
% 0. Numerical Stability Checks - Prevent state variables from going out of reasonable bounds.
if Y(1) > 50 || Y(1) < -100, warning('Vm out of range: %.2f mV', Y(1)); Y(1) = max(min(Y(1), 50), -100); end
if Y(6) < 1 || Y(6) > 200, warning('Ki out of range: %.2f mM', Y(6)); Y(6) = max(min(Y(6), 200), 1); end
if Y(7) < 0.1 || Y(7) > 150, warning('Nai out of range: %.2f mM', Y(7)); Y(7) = max(min(Y(7), 150), 0.1); end
if Y(9) < 0.001 || Y(9) > 1000, warning('Caos out of range: %.2f uM', Y(9)); Y(9) = max(min(Y(9), 1000), 0.001); end
if Y(10) < 0.001 || Y(10) > 100, warning('Cais out of range: %.2f uM', Y(10)); Y(10) = max(min(Y(10), 100), 0.001); end
if Y(18) < 0 || Y(18) > 50, Y(18) = max(min(Y(18), 50), 0); end


%% Handle the light stimulus (Smoother Gaussian Waveform with Reduced Amplitude)
% This section implements a smooth, symmetrical Gaussian stimulus.

% 1. Define Gaussian model parameters
t_peak = 1.75;         % Time of the peak stimulus (s), serving as the center of symmetry.
tau_decay_smoother = 0.6; % Time constant for a smoother curve.
sigma = tau_decay_smoother / sqrt(2);

% 2. Calculate the Gaussian waveform based on time 't'.
%    This formula generates a bell curve centered at t_peak, with a peak value of 1.
stim_waveform = exp( - (t - t_peak)^2 / (2 * sigma^2) );

% 3. Scale the waveform to the desired stimulus intensity.
% **MODIFICATION**: The peak light intensity value has been reduced from 4630 to 100.
peak_Jhv = 210;
max_stim_model = peak_Jhv * 0.43; % The peak stimulus value the model expects
params.stimulus = stim_waveform * max_stim_model;


%% 1. Unpack State Variables
states = struct();
states.Vm = Y(1); states.p_CaL = Y(2); states.p_h = Y(3); states.p_Kv = Y(4); states.m_KCa = Y(5);
states.Ki = Y(6); states.Nai = Y(7); states.Cli = Y(8); states.Caos = Y(9); states.Cais = Y(10); states.Caif = Y(11);
states.Ca_buff_os = Y(12); states.Ca_ls = Y(13); states.Ca_hs = Y(14); states.Ca_lf = Y(15); states.Ca_hf = Y(16); states.Cab = Y(17);
states.cGMP = Y(18); states.PDE_active = Y(19); states.Arr = Y(20); states.Arr_di = Y(21); states.Arr_tetra = Y(22);
states.G_GTP = Y(23); states.Ga_GDP = Y(24); states.Ga_GTP = Y(25); states.Ga_GTP_PDE_a_Ga_GTP = Y(26);
states.Ga_GTP_a_PDE_a_Ga_GTP = Y(27); states.Gbg = Y(28); states.Gt = Y(29);
states.Ops = Y(30); states.Ops_G = Y(31); states.Ops_G_GTP = Y(32); states.Ops_Gt = Y(33);
states.PDE = Y(34); states.PDE_Ga_GTP = Y(35); states.PDE_a_Ga_GTP = Y(36);
states.R = Y(37); states.R0 = Y(38); states.R0_G = Y(39); states.R0_G_GTP = Y(40);
states.R0_Gt = Y(41); states.R0_RKpre = Y(42); states.R1 = Y(43); states.R1_Arr = Y(44); states.R1_G = Y(45); states.R1_G_GTP = Y(46);
states.R1_Gt = Y(47); states.R1_RKpost = Y(48); states.R1_RKpre = Y(49);
states.R2 = Y(50); states.R2_Arr = Y(51); states.R2_G = Y(52); states.R2_G_GTP = Y(53);
states.R2_Gt = Y(54); states.R2_RKpost = Y(55); states.R2_RKpre = Y(56);
states.R3 = Y(57); states.R3_Arr = Y(58); states.R3_G = Y(59); states.R3_G_GTP = Y(60);
states.R3_Gt = Y(61); states.R3_RKpost = Y(62); states.R3_RKpre = Y(63);
states.R4 = Y(64); states.R4_Arr = Y(65); states.R4_G = Y(66); states.R4_G_GTP = Y(67);
states.R4_Gt = Y(68); states.R4_RKpost = Y(69); states.R4_RKpre = Y(70);
states.R5 = Y(71); states.R5_Arr = Y(72); states.R5_G = Y(73); states.R5_G_GTP = Y(74);
states.R5_Gt = Y(75); states.R5_RKpost = Y(76); states.R5_RKpre = Y(77);
states.R6 = Y(78); states.R6_Arr = Y(79); states.R6_G = Y(80); states.R6_G_GTP = Y(81);
states.R6_Gt = Y(82); states.R6_RKpost = Y(83); states.R6_RKpre = Y(84);
states.RGS = Y(85); states.RGS_Ga_GTP_a_PDE_a_Ga_GTP = Y(86); states.RGS_PDE_a_Ga_GTP = Y(87);
states.RK = Y(88); states.R_Gt = Y(89);
states.RecT = Y(90); states.RecR_Ca = Y(91); states.RecR_Ca_RK = Y(92);

%% 2. Calculate All Intermediate Variables and Derivatives
[dStates_photo, E_star, ~] = phototransduction_cascade(states, params.stimulus, params);
currents = calculate_currents(states, params);
I_total_pA = currents.I_total;
I_total_A = I_total_pA * 1e-12;
dVm_dt_Vs = -I_total_A / (params.Cm * 1e-12);
dVm_dt = dVm_dt_Vs * 1000;
Vm = states.Vm;
alpha_h = 80 / (1 + exp((Vm + 78) / 14));
beta_h = 180 / (1 + exp(-(Vm + 8) / 19));
dp_h_dt = (alpha_h * (1 - states.p_h) - beta_h * states.p_h);
alpha_Kv = 0.5 + (29.5) / (1 + exp(-(Vm - 5) / 17));
beta_Kv = 0.25 + (99.75) / (1 + exp((Vm + 92) / 11));
dp_Kv_dt = (alpha_Kv * (1 - states.p_Kv) - beta_Kv * states.p_Kv);
alpha_CaL = 1000 / (1 + exp(-(Vm + 20) / 10));
beta_CaL = 1000 / (1 + exp((Vm + 55) / 10));
dp_CaL_dt = (alpha_CaL * (1 - states.p_CaL) - beta_CaL * states.p_CaL);
if abs(80 - Vm) < 1e-6, alpha_mKCa = 15 / (1/40); else, alpha_mKCa = 15 * (80 - Vm) / (exp((80 - Vm) / 40) - 1); end
beta_mKCa = 20 * exp(-Vm / 35);
dm_KCa_dt = alpha_mKCa * (1 - states.m_KCa) - beta_mKCa * states.m_KCa;
J_NKCC1 = -params.k_NKCC1 * (0.1 * (1 / (1 + exp(16 - params.Ko))) * log(states.Ki * states.Cli / params.Ko / params.Clo) + log(states.Nai * states.Cli / params.Nao / params.Clo));
J_KCC2 = -params.k_KCC2 * (0.3 * log(states.Ki * states.Cli / params.Ko / params.Clo));
dK_dt = -(currents.I_CNG_K - currents.I_NCKX + currents.I_h_K + currents.I_Kv + currents.I_CaL_K + currents.I_KCa - 2 * currents.I_NaK + currents.I_L_K) * 1e-12 / (1 * params.F * params.V_cell) * 1e3 + J_NKCC1 + J_KCC2;
dNa_dt = -(currents.I_CNG_Na + 4 * currents.I_NCKX + currents.I_h_Na + currents.I_CaL_Na + 3 * currents.I_NaK + 3 * currents.I_NCX + currents.I_L_Na) * 1e-12 / (1 * params.F * params.V_cell) * 1e3 + J_NKCC1;
dCl_dt = -(currents.I_ClCa + currents.I_L_Cl) * 1e-12 / (-1 * params.F * params.V_cell) * 1e3 + 2 * J_NKCC1 + J_KCC2;
dcGMP_dt = params.alfamax / (1 + (states.Caos / params.Kc1)^params.m1) + params.alfamax / (1 + (states.Caos / params.Kc2)^params.m2) - (params.betadark + params.betasub * E_star) * states.cGMP;
J_dif = params.k3 * (states.Cais - states.Caos);
dCa_buff_os_dt = params.k1 * (params.eT - states.Ca_buff_os) * states.Caos - params.k2 * states.Ca_buff_os;
dCab_dt = dCa_buff_os_dt;
dCaos_dt = -(currents.I_CNG_Ca - 2 * currents.I_NCKX + currents.I_L_Caos) * 1e-12 / (2 * params.F * params.V_os) * 1e6 - dCa_buff_os_dt - J_dif / params.V_os;
dCa_ls_dt = params.Lb1 * states.Cais * (params.BL - states.Ca_ls) - params.Lb2 * states.Ca_ls;
dCa_hs_dt = params.Hb1 * states.Cais * (params.BH - states.Ca_hs) - params.Hb2 * states.Ca_hs;
dCa_lf_dt = params.Lb1 * states.Caif * (params.BL - states.Ca_lf) - params.Lb2 * states.Ca_lf;
dCa_hf_dt = params.Hb1 * states.Caif * (params.BH - states.Ca_hf) - params.Hb2 * states.Ca_hf;
dif = params.DCa * params.S_1 / params.delta * (states.Cais - states.Caif);
dCaif_dt = dif / params.Vif - dCa_lf_dt - dCa_hf_dt;
dCais_dt = -(currents.I_CaL_Ca + currents.I_PMCA - 2 * currents.I_NCX + currents.I_L_Cais) / (2 * params.F * params.Vis) * 1e-6 - dCa_ls_dt - dCa_hs_dt + (dif / params.Vis) + (J_dif / params.Vis);

%% 3. Pack Derivatives
dYdt = zeros(94, 1);
dYdt(1) = dVm_dt; dYdt(2) = dp_CaL_dt; dYdt(3) = dp_h_dt; dYdt(4) = dp_Kv_dt; dYdt(5) = dm_KCa_dt;
dYdt(6) = dK_dt; dYdt(7) = dNa_dt; dYdt(8) = dCl_dt; dYdt(9) = dCaos_dt; dYdt(10) = dCais_dt; dYdt(11) = dCaif_dt;
dYdt(12) = dCa_buff_os_dt; dYdt(13) = dCa_ls_dt; dYdt(14) = dCa_hs_dt; dYdt(15) = dCa_lf_dt; dYdt(16) = dCa_hf_dt; dYdt(17) = dCab_dt;
dYdt(18) = dcGMP_dt;
dYdt(19) = 0; % Legacy variable
dYdt(20:92) = cell2mat(struct2cell(dStates_photo));
dYdt(93) = 0; % Rtot is a constant parameter.
dYdt(94) = 0; % E_star is an intermediate calculated value, not a state.

end