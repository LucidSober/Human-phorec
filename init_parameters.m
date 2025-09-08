% init_parameters.m
function params = init_parameters()
    % Returns a struct containing all model parameters.
    
    %% Phototransduction Cascade Parameters
    % RK (rhodopsin kinase) Parameters
    params.kRK1_0 = 0.1724; params.omega = 2.5; params.kRK2 = 250;
    params.kRK3_ATP = 4000; params.kRK4 = 250;
    
    % Arrestin Parameters
    params.kArr = 0.0000099147; params.kA2 = 0.026; params.m_Arr = 0.0000095475;
    params.kA3 = 1.1651; params.kA4 = 0.00000029965; params.kA5 = 0.424;
    
    % Opsin Parameters
    params.kOps = 0.00000000000061172; params.kRrecyc = 0.0007;
    
    % G-protein (Transducin) Parameters
    params.omega_G = 0.6; params.kG1_0 = 0.001; params.kG2 = 2200; params.kG3 = 8500;
    params.kG4_GDP = 400; params.kG5_GTP = 3500; params.kG6 = 8500; params.kG7 = 200;
    params.kGrecyc = 2; params.kGshutoff = 0.05;
    
    % PDE Parameters
    params.kP1 = 0.05497; params.kP1_rev = 0; params.kP2 = 940.7;
    params.kP3 = 0.0000000014983; params.kP4 = 21.088; params.kPDEshutoff = 0.1;
    
    % RGS Parameters
    params.kRGS1 = 0.000048182; params.kRGS2 = 98;
    
    % Recoverin Parameters
    params.kRec1 = 0.011; params.kRec2 = 0.05; params.kRec3 = 0.00041081; params.kRec4 = 0.610084;
    
    % Thermal Activation
    params.ktherm = 0.0238; % Rate of spontaneous rhodopsin activation
    
    %% cGMP and Calcium Parameters
    params.Vcyto = 0.03916; % Cytoplasmic volume (pL)
    params.Kc1 = 0.171; % Hill coefficient for Ca binding (uM)
    params.Kc2 = 0.059; % Hill coefficient for Ca binding (uM)
    params.m1 = 3;      % Hill exponent
    params.m2 = 1.5;    % Hill exponent
    params.alfamax = 60; % Max rate of cGMP synthesis (uM/s)
    params.betadark = 3.1872748181298; % Basal cGMP hydrolysis rate (s^-1)
    params.betasub = 0.021826; % Light-driven cGMP hydrolysis rate (s^-1/molecule)
    
    % Calcium-related Parameters
    params.fCa = 0.12;           % Fraction of CNG current carried by Ca2+
    params.Jdark = 14.87;        % Dark current (pA)
    params.cGMPdark = 6.4944;    % Dark cGMP concentration (uM)
    params.ncg = 3.8;            % Hill coefficient for cGMP gating of CNG channels
    params.gammaCa = 981.3558;   % Calcium extrusion rate (s^-1)
    params.Caos_dark = 0.25;     % Outer segment Ca2+ concentration in dark (uM)
    params.Caos_0 = 0.023;       % Minimum outer segment Ca2+ concentration (uM)
    
    % Calcium Buffering Parameters
    params.k1 = 9.37059; % Forward rate for Ca buffering (uM^-1 s^-1)
    params.k2 = 46.412;  % Backward rate for Ca buffering (s^-1)
    params.eT = 400;     % Total buffer concentration (uM)
    params.k3 = 1e-20;   % Calcium diffusion rate (s^-1), likely very small
    
    %% Cell Parameters
    params.Cm = 3.6;        % Membrane capacitance (pF)
    params.F = 96485.34;    % Faraday's constant (C/mol)
    params.R_gas = 8.3145;  % Ideal gas constant (J/(mol*K))
    params.T = 310;         % Absolute temperature (K)
    params.NA = 6.022e23;   % Avogadro's number (molecules/mol)
    
    % Volume Parameters (in L)
    params.V_os = 3.916e-14; % Outer segment volume
    params.V_is = 2.3496e-14; % Inner segment volume
    params.V_if = 1.5664e-14; % Inner segment fiber volume
    params.V_cell = params.V_os + params.V_is + params.V_if;
    params.Vos = params.V_os; % Aliases for convenience
    params.Vis = params.V_is;
    params.Vif = params.V_if;
    
    %% Ion Concentrations
    % External Concentrations (mM)
    params.Cao = 1.6;  % External Ca2+
    params.Nao = 145;  % External Na+
    params.Ko = 5;     % External K+
    params.Clo = 110;  % External Cl-
    
    % Initial Internal Concentrations (for reference)
    params.Nai_0 = 5;    % Internal Na+ (mM)
    params.Ki_0 = 140;   % Internal K+ (mM)
    params.Cli_0 = 30;   % Internal Cl- (mM)
    params.Cais_0 = 0.1; % Internal Ca2+ in inner segment (uM)
    
    %% Calcium System Parameters
    params.DCa = 6e-8;      % Ca2+ diffusion coefficient (dm^2/s)
    params.delta = 3e-5;    % Diffusion distance (dm)
    params.S_1 = 3.142e-8;  % Surface area (dm^2)
    
    % Calcium Buffering Pools (Inner Segment)
    params.Lb1 = 0.4;   % Low-affinity buffer forward rate (s^-1 uM^-1)
    params.Lb2 = 0.21;  % Low-affinity buffer backward rate (s^-1)
    params.Hb1 = 100;   % High-affinity buffer forward rate (s^-1 uM^-1)
    params.Hb2 = 90;    % High-affinity buffer backward rate (s^-1)
    params.BL = 500;    % Low-affinity buffer total concentration (uM)
    params.BH = 300;    % High-affinity buffer total concentration (uM)
    
    %% Channel Conductances (nS)
    params.gKv = 0.2;     % Voltage-gated K+ conductance
    params.gCl = 2;       % Ca2+-activated Cl- conductance
    params.gKCa = 0.113;  % Ca2+-activated K+ conductance

    % Ih Channel Permeabilities (from supplementary material S3)
    params.P_K_h = 0.15;
    params.P_Na_h = 0.05;
    
    % ICaL Channel Permeabilities (from supplementary material S3)
    k_CaL_scale_factor = 0.65;
    params.P_Ca_CaL = 1 * k_CaL_scale_factor;
    params.P_K_CaL = 0.000365 * k_CaL_scale_factor;
    params.P_Na_CaL = 0.0000185 * k_CaL_scale_factor;
    
    % CNG Channel Permeabilities (from supplementary material S3)
    params.P_K_CNG = 1 / 3.6;
    params.P_Na_CNG = 1 / 3.6;
    params.P_Ca_CNG = 6 / 3.6;
    
    %% Pump and Exchanger Parameters
    params.kPMCA = 10;      % PMCA pump rate (pA/uM)
    params.Km_Cai = 0.0005; % PMCA calcium affinity (mM)
    params.kNCX = 250;      % NCX rate (pA)
    params.k_NKCC1 = 0.005; % NKCC1 cotransporter rate (s^-1)
    params.k_KCC2 = 0.046;  % KCC2 cotransporter rate (s^-1)
    params.kNKCC1 = params.k_NKCC1; % Alias
    params.kKCC2 = params.k_KCC2;   % Alias
    
    %% Leak Conductances (nS)
    params.Gleak_K = 0.00011488937611594242;
    params.Gleak_Na = 0.000071150566012260358;
    params.Gleak_Cais = 0.0002481452439628276;
    params.Gleak_Cl = 0.00035986965484864007;
    params.Gleak_Caos = 0.000021232320024139139;
    
    %% Numerical Calculation Parameters
    params.numConcFactor = 1 / (6.022e5 * params.Vcyto); % Concentration conversion factor
    
    %% Total Amounts for Phototransduction (in molecules)
    params.Rtot = 100000000.0;
    params.PDEtot = 2000000.0;
    params.Gtot = 10000000.0;
    params.RGStot = 100000.0;
    params.ArrTot = 7074600.0;
end