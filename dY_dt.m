function Y_prime = dY_dt(t, Y)

global params states_r

states_r=Y2States(Y);

%% check
if states_r.Vm > 50 || states_r.Vm < -100
    warning('Vm out of range: %.2f mV', Y(1));
    Y(1) = max(min(Y(1), 50), -100);
end

% 离子浓度限制 (mM)
if states_r.Ki < 1 || states_r.Ki > 200
    warning('Ki out of range: %.2f mM', Y(6));
    Y(6) = max(min(Y(6), 200), 1);
end

if states_r.Nai < 0.1 || states_r.Nai > 150 
    warning('Nai out of range: %.2f mM', Y(7));
    Y(7) = max(min(Y(7), 150), 0.1);
end

% 钙浓度限制 (μM)
if states_r.Caos < 0.001 || states_r.Caos > 1000
    warning('Caos out of range: %.2f μM', Y(9));
    Y(9) = max(min(Y(9), 1000), 0.001);
end
if states_r.Cais < 0.001 || states_r.Cais > 100 
    warning('Cais out of range: %.2f μM', Y(10));
    Y(10) = max(min(Y(10), 100), 0.001);
end

% cGMP限制 (μM)
if states_r.cGMP < 0 || states_r.cGMP > 50
    Y(18) = max(min(Y(18), 50), 0);
end

states_r=Y2States(Y); % update states record

% 处理光刺激
if t >= params.flash_start && t < (params.flash_start + params.flash_duration)
    params.stimulus = params.Jhv * 0.43;  % 0.43是收集面积
else
    params.stimulus = 0;
end

%% 2. 计算所有中间变量和导数
% --- 光转导级联 ---
[dStates_photo, E_star, ~] = cal_phototransduction_cascade();


% --- 离子电流 ---
ion_currents = cal_ion_currents();


% --- 膜电位导数 (dVm/dt) ---
I_total_A = ion_currents.I_total * 1e-12; % pA -> A
dVm_dt_V = -I_total_A / (params.Cm * 1e-12); % A / F -> V/s
dVm_dt = dVm_dt_V * 1000; % V/s -> mV/s


% --- 门控变量导数 (dp/dt, dm/dt) ---
Vm = states_r.Vm;
% Ih gating (使用补充材料中的正确方程)
alpha_h = 80 / (1 + exp((Vm + 78) / 14));
beta_h = 180 / (1 + exp(-(Vm + 8) / 19));
dp_h_dt = (alpha_h * (1 - states_r.p_h) - beta_h * states_r.p_h) ;

% IKv gating
alpha_Kv = 0.5 + (29.5) / (1 + exp(-(Vm - 5) / 17));
beta_Kv = 0.25 + (99.75) / (1 + exp((Vm + 92) / 11));
dp_Kv_dt = (alpha_Kv * (1 - states_r.p_Kv) - beta_Kv * states_r.p_Kv) ;

% ICaL gating
alpha_CaL = 1000 / (1 + exp(-(Vm + 20) / 10));
beta_CaL = 1000 / (1 + exp((Vm + 55) / 10));
dp_CaL_dt = (alpha_CaL * (1 - states_r.p_CaL) - beta_CaL * states_r.p_CaL) ;


% 添加数值稳定性
if abs(Vm - 5) < 1e-6
    alpha_Kv = 0.5 + (30 - 0.5) / 2;  % 使用极限值
else
    alpha_Kv = 0.5 + (30 - 0.5) / (1 + exp(-(Vm - 5) / 17));
end

if abs(Vm + 92) < 1e-6
    beta_Kv = 0.25 + (100 - 0.25) / 2;  % 使用极限值
else
    beta_Kv = 0.25 + (100 - 0.25) / (1 + exp((Vm + 92) / 11));
end

% IKCa gating - 数值稳定性检查
if abs(80 - Vm) < 1e-6
    alpha_mKCa = 15 / (1/40);  % 洛必达法则
else
    alpha_mKCa = 15 * (80 - Vm) / (exp((80 - Vm) / 40) - 1);
end
beta_mKCa = 20 * exp(-Vm / 35);

% 防止除零错误
eps_safe = 1e-9;
if abs(states_r.Ki) < eps_safe
    states_r.Ki = eps_safe;
end
if abs(states_r.Nai) < eps_safe
    states_r.Nai = eps_safe;
end
if abs(states_r.Cli) < eps_safe
    states_r.Cli = eps_safe;
end

dm_KCa_dt = alpha_mKCa * (1 - states_r.m_KCa) - beta_mKCa * states_r.m_KCa;


% --- 离子浓度和钙缓冲导数 (d[Ion]/dt) ---
J_NKCC1 = -params.k_NKCC1 * (0.1 * (1 / (1 + exp(16 - params.Ko))) * log(states_r.Ki * states_r.Cli / params.Ko / params.Clo) + ...
            log(states_r.Nai * states_r.Cli / params.Nao / params.Clo));
J_KCC2 = -params.k_KCC2 * (0.3 * log(states_r.Ki * states_r.Cli / params.Ko / params.Clo));

% 在derivatives.m中添加JNKCC1和JKCC2对总体离子浓度变化的贡献
dK_dt = -(ion_currents.I_CNG_K - ion_currents.I_NCKX + ion_currents.I_h_K + ion_currents.I_Kv + ...
      ion_currents.I_CaL_K + ion_currents.I_KCa - 2 * ion_currents.I_NaK + ion_currents.I_L_K) * 1e-12 / ...
      (1 * params.F * params.V_cell) * 1e3 + J_NKCC1 + J_KCC2; % 添加了共转运体的贡献

dNa_dt = -(ion_currents.I_CNG_Na + 4 * ion_currents.I_NCKX + ion_currents.I_h_Na + ...
       ion_currents.I_CaL_Na + 3 * ion_currents.I_NaK + 3 * ion_currents.I_NCX + ion_currents.I_L_Na) * 1e-12 / ...
       (1 * params.F * params.V_cell) * 1e3 + J_NKCC1; % 添加了NKCC1的贡献

dCl_dt = -(ion_currents.I_ClCa + ion_currents.I_L_Cl) * 1e-12 / ...
       (-1 * params.F * params.V_cell) * 1e3 + 2 * J_NKCC1 + J_KCC2; % 已包含

% cGMP 导数
dcGMP_dt = params.alfamax / (1 + (states_r.Caos / params.Kc1)^params.m1) + ...
           params.alfamax / (1 + (states_r.Caos / params.Kc2)^params.m2) - ...
           (params.betadark + params.betasub * E_star) * states_r.cGMP;

% 钙离子导数 (包含缓冲和扩散)
J_dif = params.k3 * (states_r.Cais - states_r.Caos);
dCa_buff_os_dt = params.k1 * (params.eT - states_r.Ca_buff_os) * states_r.Caos - params.k2 * states_r.Ca_buff_os;
dCab_dt = dCa_buff_os_dt; 
dCaos_dt = -(ion_currents.I_CNG_Ca - 2 * ion_currents.I_NCKX + ion_currents.I_L_Caos) * 1e-12 / ...
            (2 * params.F * params.V_os) * 1e6 - dCa_buff_os_dt - J_dif / params.V_os; % mol/s -> uM/s

% 内段钙和其缓冲 (修正了所有参数名以匹配 init_parameters.m)
dCa_ls_dt = params.Lb1 * states_r.Cais * (params.BL - states_r.Ca_ls) - params.Lb2 * states_r.Ca_ls;
dCa_hs_dt = params.Hb1 * states_r.Cais * (params.BH - states_r.Ca_hs) - params.Hb2 * states_r.Ca_hs;
dCa_lf_dt = params.Lb1 * states_r.Caif * (params.BL - states_r.Ca_lf) - params.Lb2 * states_r.Ca_lf;
dCa_hf_dt = params.Hb1 * states_r.Caif * (params.BH - states_r.Ca_hf) - params.Hb2 * states_r.Ca_hf;
dif = params.DCa * params.S_1 / params.delta * (states_r.Cais - states_r.Caif);
dCaif_dt = dif / params.Vif - dCa_lf_dt - dCa_hf_dt;
dCais_dt = -(ion_currents.I_CaL_Ca + ion_currents.I_PMCA - 2 * ion_currents.I_NCX + ion_currents.I_L_Cais) / ...
(2 * params.F * params.Vis) * 1e-6 - dCa_ls_dt - dCa_hs_dt + (dif / params.Vis) + (J_dif / params.Vis);


%% 3. 打包导数 (Pack Derivatives)
% 按照与解包时完全相同的顺序，将所有导数打包到 dYdt 向量中
Y_prime = zeros(92, 1);

% 电生理状态
Y_prime(1) = dVm_dt; 
Y_prime(2) = dp_CaL_dt; 
Y_prime(3) = dp_h_dt; 
Y_prime(4) = dp_Kv_dt; 
Y_prime(5) = dm_KCa_dt;

% 离子浓度
Y_prime(6) = dK_dt; 
Y_prime(7) = dNa_dt; 
Y_prime(8) = dCl_dt; 
Y_prime(9) = dCaos_dt; 
Y_prime(10) = dCais_dt; 
Y_prime(11) = dCaif_dt;

% 钙缓冲
Y_prime(12) = dCa_buff_os_dt; 
Y_prime(13) = dCa_ls_dt; 
Y_prime(14) = dCa_hs_dt; 
Y_prime(15) = dCa_lf_dt; 
Y_prime(16) = dCa_hf_dt; 
Y_prime(17) = dCab_dt;

% 光转导状态
Y_prime(18) = dcGMP_dt;
Y_prime(19) = 0; % 遗留变量
Y_prime(20:92) = cell2mat(struct2cell(dStates_photo)); % 快速打包光转导状态导

Y_prime(93) = 0; % Rtot is calculated
Y_prime(94) = 0; % E_star is calculated




end

%[appendix]{"version":"1.0"}
%---
