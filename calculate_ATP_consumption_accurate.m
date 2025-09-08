% calculate_ATP_consumption_accurate.m
function [ATP_total, ATP_Na, ATP_Ca, ATP_photo] = calculate_ATP_consumption_accurate(currents, states, params, stimulus)
% Accurately calculates ATP consumption from different cellular processes.
% Inputs:
%   currents - A struct containing all current components.
%   states   - A struct containing all state variables.
%   params   - A struct containing all model parameters.
%   stimulus - The current light stimulus intensity.
% Outputs:
%   ATP_total - Total ATP consumption rate (molecules/s).
%   ATP_Na    - ATP consumption rate by the Na+/K+ pump (molecules/s).
%   ATP_Ca    - ATP consumption rate by the PMCA pump (molecules/s).
%   ATP_photo - ATP consumption rate by the phototransduction cascade (molecules/s).

    % Extract constants
    F = params.F;
    Vcell = params.V_cell;
    Vis = params.V_is;
    NA = params.NA;
    
    %% 1. ATP for Na+ Pumping (via Na+/K+ Pump)
    % The Na+/K+ pump consumes 1 ATP to pump out 3 Na+. Its rate matches the total Na+ influx.
    Na_influx_currents_pA = currents.I_CNG_Na + currents.I_h_Na + currents.I_CaL_Na + currents.I_L_Na;
    % Note: I_NCKX and I_NCX contribute to Na+ flux, but their stoichiometry is complex.
    % A more direct method is to use the Na/K pump current itself.
    % 1 ATP per 3 Na+ ions pumped out, which corresponds to a net charge of +1e moved out.
    % The pump current I_NaK is defined as positive outward current.
    % I_NaK = I_pump * (3 * q_e - 2 * q_e) = I_pump * q_e, where I_pump is pump cycles/sec.
    % Number of Na+ ions/sec = 3 * I_pump = 3 * (I_NaK_pA * 1e-12) / (1.602e-19 C)
    % ATP molecules/sec = (Number of Na+ ions/sec) / 3 = (I_NaK_pA * 1e-12) / (1.602e-19 C)
    
    I_NaK_Amps = currents.I_NaK * 1e-12;
    charge_per_ion = 1.60217663e-19; % Elementary charge in Coulombs
    ATP_Na = I_NaK_Amps / charge_per_ion; % This assumes 1 net charge per ATP.
    
    %% 2. ATP for Ca2+ Pumping (via PMCA)
    % The PMCA pump consumes 1 ATP to pump out 1 Ca2+.
    % The current is defined as positive outward.
    I_PMCA_Amps = currents.I_PMCA * 1e-12;
    charge_per_Ca = 2 * charge_per_ion;
    ATP_Ca = I_PMCA_Amps / charge_per_Ca;
    
    %% 3. ATP Consumption by the Phototransduction Cascade
    ATP_photo = 0; % Default to zero in darkness
    if stimulus > 0
        % In the presence of light, ATP is consumed for rhodopsin phosphorylation and G-protein activation.
        % NOTE: This section re-calculates reaction rates, which is a minor redundancy.
        % A more streamlined approach would be to pass the 'v' struct from the main ODE function.
        v_r3_0 = params.kRK3_ATP * states.R0_RKpre;
        v_r3_1 = params.kRK3_ATP * states.R1_RKpre;
        v_r3_2 = params.kRK3_ATP * states.R2_RKpre;
        v_r3_3 = params.kRK3_ATP * states.R3_RKpre;
        v_r3_4 = params.kRK3_ATP * states.R4_RKpre;
        v_r3_5 = params.kRK3_ATP * states.R5_RKpre;
        v_r10 = params.kG5_GTP * states.Ops_G;
        v_r15_0 = params.kG5_GTP * states.R0_G;
        v.r15_1 = params.kG5_GTP * states.R1_G; v.r15_2 = params.kG5_GTP * states.R2_G;
        v.r15_3 = params.kG5_GTP * states.R3_G; v.r15_4 = params.kG5_GTP * states.R4_G;
        v.r15_5 = params.kG5_GTP * states.R5_G; v.r15_6 = params.kG5_GTP * states.R6_G;
        
        % Sum of all ATP-consuming reaction rates (in molecules/sec)
        ATP_photo = v_r3_0 + v_r3_1 + v_r3_2 + v_r3_3 + v_r3_4 + v_r3_5 + ...
                    v_r10 + v_r15_0 + v.r15_1 + v.r15_2 + v.r15_3 + v.r15_4 + v.r15_5 + v.r15_6;
    end
    
    %% 4. Total ATP Consumption
    ATP_total = ATP_Na + ATP_Ca + ATP_photo;
end