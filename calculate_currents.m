% calculate_currents.m (Final Corrected Version)
function currents = calculate_currents(states, params, Vm)
% This version ensures it uses the true transmembrane potential Vm passed as an argument.

%% 1. Extract Variables
    F = params.F;
    R = params.R_gas;
    T = params.T;
    Cm = params.Cm;
    
    % CRITICAL: The line "Vm = states.Vm;" has been REMOVED.
    % The function now uses the Vm that is passed in as the third argument.
    
    Ki = states.Ki;
    Nai = states.Nai;
    Cli = states.Cli;
    Caos = states.Caos;
    Cais = states.Cais;
    cGMP = states.cGMP;
    
    p_h = states.p_h;
    p_Kv = states.p_Kv;
    p_CaL = states.p_CaL;
    m_KCa = states.m_KCa;
    
    Ko = params.Ko;
    Nao = params.Nao;
    Cao_mM = params.Cao;
    Clo = params.Clo;
    
    eps = 1e-9;

%% 2. Calculate Nernst Equilibrium Potentials (in mV)
    E_K = (R * T / F) * log(Ko / (Ki + eps)) * 1e3;
    E_Na = (R * T / F) * log(Nao / (Nai + eps)) * 1e3;
    E_Cl = (R * T / (-1 * F)) * log(Clo / (Cli + eps)) * 1e3;
    E_Caos = (R * T / (2 * F)) * log(Cao_mM / (Caos*1e-3 + eps)) * 1e3;
    E_Cais = (R * T / (2 * F)) * log(Cao_mM / (Cais*1e-3 + eps)) * 1e3;

%% 3. Calculate GHK Driving Force Terms
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
    
%% 4. Calculate Individual Ionic Currents (in pA)
    pOpen_CNG = -2 / (2 + params.fCa) * params.Jdark * (cGMP / params.cGMPdark)^params.ncg / ...
                (params.P_K_CNG * CF_K + params.P_Na_CNG * CF_Na + params.P_Ca_CNG * CF_Caos) / Cm;
    currents.I_CNG_K = params.P_K_CNG * CF_K * pOpen_CNG * Cm;
    currents.I_CNG_Na = params.P_Na_CNG * CF_Na * pOpen_CNG * Cm;
    currents.I_CNG_Ca = params.P_Ca_CNG * CF_Caos * pOpen_CNG * Cm;
    currents.I_CNG = currents.I_CNG_K + currents.I_CNG_Na + currents.I_CNG_Ca;
    Ca2_frac = (Caos - params.Caos_0) / (params.Caos_dark - params.Caos_0);
    currents.I_NCKX = -params.fCa / (params.fCa + 2) * Ca2_frac * params.Jdark;
    currents.I_h_K = params.P_K_h * CF_K * p_h * Cm;
    currents.I_h_Na = params.P_Na_h * CF_Na * p_h * Cm;
    currents.I_h = currents.I_h_K + currents.I_h_Na;
    currents.I_Kv = params.gKv * p_Kv * (Vm - E_K);
    currents.I_CaL_Ca = params.P_Ca_CaL * CF_Cais * p_CaL * Cm;
    currents.I_CaL_K = params.P_K_CaL * CF_K * p_CaL * Cm;
    currents.I_CaL_Na = params.P_Na_CaL * CF_Na * p_CaL * Cm;
    currents.I_CaL = currents.I_CaL_Ca + currents.I_CaL_K + currents.I_CaL_Na;
    mCl = 1 / (1 + exp((0.37 - Cais) / 0.09));
    currents.I_ClCa = params.gCl * mCl * (Vm - E_Cl);
    mKCas = Cais / (Cais + 0.3);
    currents.I_KCa = params.gKCa * m_KCa^2 * mKCas * (Vm - E_K);
    currents.I_PMCA = params.kPMCA * Cm * Cais / (Cais + params.Km_Cai * 1000);
    eta_NCX = 0.35; ksat_NCX = 0.1; Km_Na = 87.5; Km_Ca = 1.38;
    fNCX_Nao = 1 / (Km_Na^3 + Nao^3);
    fNCX_Cao = 1 / (Km_Ca + Cao_mM);
    fNCX_Vm = 1 / (1 + ksat_NCX * exp((eta_NCX - 1) * F * Vm_V / (R*T))) * ...
              (exp(eta_NCX * F * Vm_V / (R*T)) * Nai^3 * Cao_mM - ...
               exp((eta_NCX - 1) * F * Vm_V / (R*T)) * Nao^3 * Cais * 1e-3);
    currents.I_NCX = Cm * params.kNCX * fNCX_Nao * fNCX_Cao * fNCX_Vm;
    currents.I_NaK = 10.34 * (0.33 + (1 - 0.33) / (1 + exp(-(Vm + 57.1) / 27.6)));
    currents.I_L_K = params.Gleak_K * (Vm - E_K);
    currents.I_L_Na = params.Gleak_Na * (Vm - E_Na);
    currents.I_L_Cais = params.Gleak_Cais * (Vm - E_Cais);
    currents.I_L_Cl = params.Gleak_Cl * (Vm - E_Cl);
    currents.I_L_Caos = params.Gleak_Caos * (Vm - E_Caos);
    
%% 5. Total Current
    currents.I_total = currents.I_CNG + currents.I_NCKX + currents.I_h + currents.I_Kv + ...
                       currents.I_CaL + currents.I_ClCa + currents.I_KCa + currents.I_NCX + ...
                       currents.I_PMCA + currents.I_NaK + currents.I_L_K + currents.I_L_Na + ...
                       currents.I_L_Cais + currents.I_L_Cl + currents.I_L_Caos;
end