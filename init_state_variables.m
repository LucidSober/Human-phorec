% init_state_variables.m
function states = init_state_variables()
    % Initializes all state variables into a unified struct.
    
    %% Membrane Potential and Gating Variables
    states.Vm = -36.186; % mV
    states.p_h   = 0.103; % Gating variable for I_h
    states.p_Kv  = 0.769; % Gating variable for I_Kv
    states.p_CaL = 0.556; % Gating variable for I_CaL
    states.m_KCa = 0.642; % Gating variable for I_KCa
    
    %% Ion Concentrations (mM)
    states.Ki = 140; 
    states.Nai = 5;  
    states.Cli = 30; 
    
    %% Calcium System (uM)
    states.cGMP = 6.5;         % cGMP concentration
    states.PDE_active = 0;     % Legacy state variable
    states.Caos = 0.25;        % Outer segment free Ca2+
    states.Cais = 0.1;         % Inner segment free Ca2+
    states.Caif = 0.1;         % Inner segment fiber free Ca2+
    states.Cab = 19.21989;     % Outer segment buffered Ca2+
    
    % Inner Segment Calcium Buffers (names corrected for consistency)
    states.Ca_ls = 80; % Ca2+ bound to low-affinity sites in inner segment (uM)
    states.Ca_hs = 30; % Ca2+ bound to high-affinity sites in inner segment (uM)
    states.Ca_lf = 80; % Ca2+ bound to low-affinity sites in fiber (uM)
    states.Ca_hf = 30; % Ca2+ bound to high-affinity sites in fiber (uM)
    
    %% Phototransduction Cascade State Variables (molecules)
    % These values represent the dark-adapted steady state.
    states.Arr = 1260752.939;
    states.Arr_di = 1123332.707;
    states.Arr_tetra = 891795.411;
    states.G_GTP = 0;
    states.Ga_GDP = 0;
    states.Ga_GTP = 0;
    states.Ga_GTP_PDE_a_Ga_GTP = 0;
    states.Ga_GTP_a_PDE_a_Ga_GTP = 0;
    states.Gbg = 0;
    states.Gt = 8152518.874;
    states.Ops = 0;
    states.Ops_G = 0;
    states.Ops_G_GTP = 0;
    states.Ops_Gt = 0;
    states.PDE = 2000000.0;
    states.PDE_Ga_GTP = 0;
    states.PDE_a_Ga_GTP = 0;
    states.R = 98152518.874;
    states.R_Gt = 1847481.125;
    states.RGS = 100000.0;
    states.RK = 579.637;
    states.RecR_Ca = 510930.691;
    states.RecR_Ca_RK = 199420.362;
    states.RecT = 9289648.942;
    
    % R0-R6 states and their associated complexes (initially zero for a dark-adapted state)
    states.R0 = 0; states.R0_G = 0; states.R0_G_GTP = 0; states.R0_Gt = 0; states.R0_RKpre = 0;
    states.R1 = 0; states.R1_Arr = 0; states.R1_G = 0; states.R1_G_GTP = 0; 
    states.R1_Gt = 0; states.R1_RKpost = 0; states.R1_RKpre = 0;
    states.R2 = 0; states.R2_Arr = 0; states.R2_G = 0; states.R2_G_GTP = 0; 
    states.R2_Gt = 0; states.R2_RKpost = 0; states.R2_RKpre = 0;
    states.R3 = 0; states.R3_Arr = 0; states.R3_G = 0; states.R3_G_GTP = 0; 
    states.R3_Gt = 0; states.R3_RKpost = 0; states.R3_RKpre = 0;
    states.R4 = 0; states.R4_Arr = 0; states.R4_G = 0; states.R4_G_GTP = 0; 
    states.R4_Gt = 0; states.R4_RKpost = 0; states.R4_RKpre = 0;
    states.R5 = 0; states.R5_Arr = 0; states.R5_G = 0; states.R5_G_GTP = 0; 
    states.R5_Gt = 0; states.R5_RKpost = 0; states.R5_RKpre = 0;
    states.R6 = 0; states.R6_Arr = 0; states.R6_G = 0; states.R6_G_GTP = 0; 
    states.R6_Gt = 0; states.R6_RKpost = 0; states.R6_RKpre = 0;
    
    states.RGS_Ga_GTP_a_PDE_a_Ga_GTP = 0;
    states.RGS_PDE_a_Ga_GTP = 0;
    
    % Calculated value, not a state, but initialized here
    states.E_star = 0;

    % Parameter, but initialized here for struct completeness
    states.Rtot = 100000000.0;
end