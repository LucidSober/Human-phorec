%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function ion_currents = cal_ion_currents()
% CALCULATE_CURRENTS 计算模型中所有的离子电流
%   此版本修正了函数接口以匹配derivatives.m的调用

global params states_r

%% 1. 提取变量
    F = params.F;
    R = params.R_gas;
    T = params.T;
    Cm = params.Cm;
    
    % 从states结构体中提取所有需要的变量
    Vm = states_r.Vm;
    Ki = states_r.Ki;
    Nai = states_r.Nai;
    Cli = states_r.Cli;
    Caos = states_r.Caos;
    Cais = states_r.Cais;
    cGMP = states_r.cGMP; % 现在可以正确获取
    
    p_h = states_r.p_h;
    p_Kv = states_r.p_Kv;
    p_CaL = states_r.p_CaL;
    m_KCa = states_r.m_KCa;
    
    Ko = params.Ko;
    Nao = params.Nao;
    Cao_mM = params.Cao; % 假设params.Cao已经是mM单位
    Clo = params.Clo;
    
    eps = 1e-9;

%% 2. 计算能斯特(Nernst)平衡电位 (mV)
    E_K = (R * T / F) * log(Ko / (Ki + eps)) * 1e3;
    E_Na = (R * T / F) * log(Nao / (Nai + eps)) * 1e3;
    E_Cl = (R * T / (-1 * F)) * log(Clo / (Cli + eps)) * 1e3;
    E_Caos = (R * T / (2 * F)) * log(Cao_mM / (Caos*1e-3 + eps)) * 1e3;
    E_Cais = (R * T / (2 * F)) * log(Cao_mM / (Cais*1e-3 + eps)) * 1e3;

%% 3. 计算GHK驱动项
    Vm_V = Vm * 1e-3;
    if abs(Vm_V) < 1e-7
        CF_K = F/(R*T) * (Ki - Ko);
        CF_Na = F/(R*T) * (Nai - Nao);
        CF_Caos = 2*F/(R*T) * (Caos*1e-3 - Cao_mM);
        CF_Cais = 2*F/(R*T) * (Cais*1e-3 - Cao_mM);
    else
        zFVRT_K = 1 * F * Vm_V / (R * T);
        zFVRT_Na = 1 * F * Vm_V / (R * T);
        zFVRT_Ca = 2 * F * Vm_V / (R * T);
        CF_K = zFVRT_K * (Ki - Ko * exp(-zFVRT_K)) / (1 - exp(-zFVRT_K));
        CF_Na = zFVRT_Na * (Nai - Nao * exp(-zFVRT_Na)) / (1 - exp(-zFVRT_Na));
        CF_Caos = zFVRT_Ca * (Caos*1e-3 - Cao_mM * exp(-zFVRT_Ca)) / (1 - exp(-zFVRT_Ca));
        CF_Cais = zFVRT_Ca * (Cais*1e-3 - Cao_mM * exp(-zFVRT_Ca)) / (1 - exp(-zFVRT_Ca));
    end
    
%% 4. 计算各个离子电流 (pA)
    % --- 1. I_CNG ---
    pOpen_CNG = -2 / (2 + params.fCa) * params.Jdark * (cGMP / params.cGMPdark)^params.ncg / ...
                (params.P_K_CNG * CF_K + params.P_Na_CNG * CF_Na + params.P_Ca_CNG * CF_Caos) / Cm;
    
    ion_currents.I_CNG_K = params.P_K_CNG * CF_K * pOpen_CNG * Cm;
    ion_currents.I_CNG_Na = params.P_Na_CNG * CF_Na * pOpen_CNG * Cm;
    ion_currents.I_CNG_Ca = params.P_Ca_CNG * CF_Caos * pOpen_CNG * Cm;
    ion_currents.I_CNG = ion_currents.I_CNG_K + ion_currents.I_CNG_Na + ion_currents.I_CNG_Ca;
    
    % --- 2. I_NCKX ---
    Ca2_frac = (Caos - params.Caos_0) / (params.Caos_dark - params.Caos_0);
    ion_currents.I_NCKX = -params.fCa / (params.fCa + 2) * Ca2_frac * params.Jdark;
    
    % --- 3. I_h ---
    ion_currents.I_h_K = params.P_K_h * CF_K * p_h * Cm;
    ion_currents.I_h_Na = params.P_Na_h * CF_Na * p_h * Cm;
    ion_currents.I_h = ion_currents.I_h_K + ion_currents.I_h_Na;
    
    % --- 4. I_Kv ---
    ion_currents.I_Kv = params.gKv * p_Kv * (Vm - E_K);
    
    % --- 5. I_CaL ---
    ion_currents.I_CaL_Ca = params.P_Ca_CaL * CF_Cais * p_CaL * Cm;
    ion_currents.I_CaL_K = params.P_K_CaL * CF_K * p_CaL * Cm;
    ion_currents.I_CaL_Na = params.P_Na_CaL * CF_Na * p_CaL * Cm;
    ion_currents.I_CaL = ion_currents.I_CaL_Ca + ion_currents.I_CaL_K + ion_currents.I_CaL_Na;
    
    % --- 6. I_ClCa ---
    mCl = 1 / (1 + exp((0.37 - Cais) / 0.09));
    ion_currents.I_ClCa = params.gCl * mCl * (Vm - E_Cl);
    
    % --- 7. I_KCa ---
    mKCas = Cais / (Cais + 0.3);
    ion_currents.I_KCa = params.gKCa * m_KCa^2 * mKCas * (Vm - E_K);
    
    % --- 8. I_PMCA ---
    ion_currents.I_PMCA = params.kPMCA * Cm * Cais / (Cais + params.Km_Cai * 1000);

    % --- 9. I_NCX --- 
    eta_NCX  = 0.35; ksat_NCX = 0.1; Km_Na = 87.5; Km_Ca = 1.38;
    fNCX_Nao = 1 / (Km_Na^3 + Nao^3);
    fNCX_Cao = 1 / (Km_Ca + Cao_mM);
    fNCX_Vm  = 1 / (1 + ksat_NCX * exp((eta_NCX - 1) * F * Vm_V / (R*T))) * ...
              (exp(eta_NCX * F * Vm_V / (R*T)) * Nai^3 * Cao_mM - ...
               exp((eta_NCX - 1) * F * Vm_V / (R*T)) * Nao^3 * Cais * 1e-3);
    ion_currents.I_NCX = Cm * params.kNCX * fNCX_Nao * fNCX_Cao * fNCX_Vm;
    
    % --- 10. I_NaK --- @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ion_currents.I_NaK =   10.34 * (0.33 + (1 - 0.33) / (1 + exp(-(Vm + 57.1) / 27.6)))  ;  %@@@
    
    % --- 11. 泄漏电流 ---
    ion_currents.I_L_K = params.Gleak_K * (Vm - E_K);
    ion_currents.I_L_Na = params.Gleak_Na * (Vm - E_Na);
    ion_currents.I_L_Cais = params.Gleak_Cais * (Vm - E_Cais);
    ion_currents.I_L_Cl = params.Gleak_Cl * (Vm - E_Cl);
    ion_currents.I_L_Caos = params.Gleak_Caos * (Vm - E_Caos);
    
    %% 5. 总电流 pA
    ion_currents.I_total = ion_currents.I_CNG + ion_currents.I_NCKX + ion_currents.I_h + ion_currents.I_Kv + ...
                           ion_currents.I_CaL + ion_currents.I_ClCa + ion_currents.I_KCa + ion_currents.I_NCX + ...
                           ion_currents.I_PMCA + ion_currents.I_NaK + ion_currents.I_L_K + ion_currents.I_L_Na + ...
                           ion_currents.I_L_Cais + ion_currents.I_L_Cl + ion_currents.I_L_Caos;
end

%[appendix]{"version":"1.0"}
%---
