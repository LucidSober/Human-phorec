% main.m
% This script runs the complete rod photoreceptor model simulation.
% It uses the ode15s solver for numerical stability, which is suitable for stiff ODEs.

clear all; close all; clc;

%% 1. Initialization
fprintf('1. Initializing parameters and state variables...\n');
params = init_parameters();      
states = init_state_variables();  
Y0 = pack_states_to_vector(states); 

%% 2. Simulation Setup
% Define light intensities and add to params struct for universal access
light_intensities = [1.7, 4.8, 15.2, 39.4, 125, 444, 1406, 4630];
params.light_intensities = light_intensities;
params.flash_start = 0;       % Initial default value
params.flash_duration = 0;    % Initial default value
params.stimulus = 0;          % Initial stimulus is zero (darkness)


% ODE solver options - using stricter tolerances for accuracy
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, ...
                 'MaxStep', 0.001, ...
                 'InitialStep', 1e-6);

%% 3. Find the Steady State in Darkness
fprintf('2. Finding the steady state in dark conditions...\n');
params.stimulus = 0; % Ensure dark condition
params.Jhv = 0;

tspan_ss = [0 5];
[~, Y_ss] = ode15s(@(t,Y) derivatives(t, Y, params), tspan_ss, Y0, options);
Y_steady_state = Y_ss(end, :)';

% Unpack steady state with a consistent variable name for later use
states_steady_state = unpack_vector_to_states(Y_steady_state);

% Validate the steady state against expected values
fprintf('\nSteady-State Check:\n');
fprintf('  Vm = %.2f mV (Expected: ~-36.2 mV)\n', states_steady_state.Vm);
fprintf('  cGMP = %.2f uM (Expected: ~6.5 uM)\n', states_steady_state.cGMP);
fprintf('  Caos = %.3f uM (Expected: ~0.25 uM)\n', states_steady_state.Caos);
fprintf('  Cais = %.3f uM (Expected: ~0.1 uM)\n', states_steady_state.Cais);
fprintf('  Ki = %.1f mM (Expected: ~140 mM)\n', states_steady_state.Ki);
fprintf('  Nai = %.1f mM (Expected: ~5 mM)\n', states_steady_state.Nai);
currents_steady_state = calculate_currents(states_steady_state, params);
fprintf('  Dark Current = %.2f pA (Expected: ~-15 pA)\n', currents_steady_state.I_CNG + currents_steady_state.I_NCKX);
fprintf('   Steady state achieved!\n\n');


% =========================================================================
%% 4. PART 1: Run FLASH SIMULATIONS for FIGURE 2
% =========================================================================
fprintf('3. Starting FLASH response simulations (for Figure 2)...\n');
results_flash = cell(length(params.light_intensities), 1);
t_results_flash = cell(length(params.light_intensities), 1);

t_total_flash = 5.0; 
params.flash_start = 0.5;
params.flash_duration = 0.02;

for i = 1:length(params.light_intensities)
    J_hv = params.light_intensities(i);
    fprintf('   Simulating flash %d/%d: Intensity = %.1f\n', i, length(params.light_intensities), J_hv);
    params.Jhv = J_hv;
    [t_full, Y_full] = ode15s(@(t,Y) derivatives(t, Y, params), [0, t_total_flash], Y_steady_state, options);
    t_results_flash{i} = t_full;
    results_flash{i} = unpack_solution_to_results(t_full, Y_full, params, J_hv);
end
fprintf('   Flash simulations complete.\n\n');


% =========================================================================
%% 5. PART 2: Run CONSTANT LIGHT SIMULATIONS for FIGURES 3 & 5
% =========================================================================
fprintf('4. Starting CONSTANT LIGHT simulations (for Figures 3 & 5)...\n');
results_constant = cell(length(params.light_intensities), 1);
t_results_constant = cell(length(params.light_intensities), 1);

t_total_const = 8.0; 
params.flash_start = 1.0; 
params.flash_duration = 5.0;

% Declare variables to store raw final simulation data for Figure 3
Y_final_high_intensity = [];
t_final_high_intensity = [];

for i = 1:length(params.light_intensities)
    J_hv = params.light_intensities(i);
    fprintf('   Simulating constant light %d/%d: Intensity = %.1f\n', i, length(params.light_intensities), J_hv);
    params.Jhv = J_hv;
    [t_full, Y_full] = ode15s(@(t,Y) derivatives(t, Y, params), [0, t_total_const], Y_steady_state, options);
    t_results_constant{i} = t_full;
    results_constant{i} = unpack_solution_to_results(t_full, Y_full, params, J_hv);
    
    % Save the raw data from the highest intensity simulation for Figure 3
    if i == length(params.light_intensities)
        Y_final_high_intensity = Y_full;
        t_final_high_intensity = t_full;
    end
end
fprintf('   Constant light simulations complete.\n\n');

% Create a complete light_states struct from the raw data for Figure 3
light_end_time = params.flash_start + params.flash_duration;
light_end_idx = find(t_final_high_intensity >= light_end_time, 1);
Y_light_steady_vector = Y_final_high_intensity(light_end_idx, :)';
light_states_steady = unpack_vector_to_states(Y_light_steady_vector);


%% 6. Plotting All Figures
fprintf('5. Plotting all results...\n');
% Plot Figure 2 (Flash responses)
plot_final_results_figure2(results_flash, t_results_flash, params);

% --- FIX: Call plot_figure3 with the correct, complete state structs ---
% Plot Figure 3 (ATP Bar Chart)
plot_figure3(light_states_steady, states_steady_state, params);

% Plot Figure 5 (Constant Light Analysis)
plot_figure5(results_constant, states_steady_state, params.light_intensities);

fprintf('Simulations and plotting complete!\n');

% NOTE: The following helper functions are duplicated in derivatives.m.
% For a cleaner implementation, they should be defined only once.

function Y0 = pack_states_to_vector(states)
    % Packs the states struct into a single column vector Y0.
    % The order must be consistent with unpack_vector_to_states.
    
    Y0 = zeros(94, 1);
    
    % Electrophysiological States (Gating)
    Y0(1) = states.Vm;
    Y0(2) = states.p_CaL;
    Y0(3) = states.p_h;
    Y0(4) = states.p_Kv;
    Y0(5) = states.m_KCa;

    % Ion Concentrations
    Y0(6) = states.Ki;
    Y0(7) = states.Nai;
    Y0(8) = states.Cli;
    Y0(9) = states.Caos;
    Y0(10) = states.Cais;
    Y0(11) = states.Caif;

    % Calcium Buffering States
    Y0(12) = states.Cab; % Note: Cab and Ca_buff_os are the same state.
    Y0(13) = states.Ca_ls;
    Y0(14) = states.Ca_hs;
    Y0(15) = states.Ca_lf;
    Y0(16) = states.Ca_hf;
    Y0(17) = states.Cab; % Placeholder for consistent indexing.

    % Phototransduction States
    Y0(18) = states.cGMP;
    Y0(19) = states.PDE_active; % Legacy placeholder
    % ... (packing all 73 phototransduction states)
    Y0(20) = states.Arr; Y0(21) = states.Arr_di; Y0(22) = states.Arr_tetra;
    Y0(23) = states.G_GTP; Y0(24) = states.Ga_GDP; Y0(25) = states.Ga_GTP;
    Y0(26) = states.Ga_GTP_PDE_a_Ga_GTP; Y0(27) = states.Ga_GTP_a_PDE_a_Ga_GTP;
    Y0(28) = states.Gbg; Y0(29) = states.Gt;
    Y0(30) = states.Ops; Y0(31) = states.Ops_G; Y0(32) = states.Ops_G_GTP; Y0(33) = states.Ops_Gt;
    Y0(34) = states.PDE; Y0(35) = states.PDE_Ga_GTP; Y0(36) = states.PDE_a_Ga_GTP;
    Y0(37) = states.R; Y0(38) = states.R0; Y0(39) = states.R0_G; Y0(40) = states.R0_G_GTP;
    Y0(41) = states.R0_Gt; Y0(42) = states.R0_RKpre;
    Y0(43) = states.R1; Y0(44) = states.R1_Arr; Y0(45) = states.R1_G; Y0(46) = states.R1_G_GTP;
    Y0(47) = states.R1_Gt; Y0(48) = states.R1_RKpost; Y0(49) = states.R1_RKpre;
    Y0(50) = states.R2; Y0(51) = states.R2_Arr; Y0(52) = states.R2_G; Y0(53) = states.R2_G_GTP;
    Y0(54) = states.R2_Gt; Y0(55) = states.R2_RKpost; Y0(56) = states.R2_RKpre;
    Y0(57) = states.R3; Y0(58) = states.R3_Arr; Y0(59) = states.R3_G; Y0(60) = states.R3_G_GTP;
    Y0(61) = states.R3_Gt; Y0(62) = states.R3_RKpost; Y0(63) = states.R3_RKpre;
    Y0(64) = states.R4; Y0(65) = states.R4_Arr; Y0(66) = states.R4_G; Y0(67) = states.R4_G_GTP;
    Y0(68) = states.R4_Gt; Y0(69) = states.R4_RKpost; Y0(70) = states.R4_RKpre;
    Y0(71) = states.R5; Y0(72) = states.R5_Arr; Y0(73) = states.R5_G; Y0(74) = states.R5_G_GTP;
    Y0(75) = states.R5_Gt; Y0(76) = states.R5_RKpost; Y0(77) = states.R5_RKpre;
    Y0(78) = states.R6; Y0(79) = states.R6_Arr; Y0(80) = states.R6_G; Y0(81) = states.R6_G_GTP;
    Y0(82) = states.R6_Gt; Y0(83) = states.R6_RKpost; Y0(84) = states.R6_RKpre;
    Y0(85) = states.RGS; Y0(86) = states.RGS_Ga_GTP_a_PDE_a_Ga_GTP; Y0(87) = states.RGS_PDE_a_Ga_GTP;
    Y0(88) = states.RK; Y0(89) = states.R_Gt;
    Y0(90) = states.RecT; Y0(91) = states.RecR_Ca; Y0(92) = states.RecR_Ca_RK;
    Y0(93) = states.Rtot;   % This is a parameter, but packed for consistency.
    Y0(94) = states.E_star; % This is a calculated value, but packed for consistency.
end

function states = unpack_vector_to_states(Y)
    % Unpacks a single time point column vector Y back into the states struct.
    % The order must be consistent with pack_states_to_vector.
    
    states = struct();
    
    % Electrophysiological States (Gating)
    states.Vm = Y(1); states.p_CaL = Y(2); states.p_h = Y(3); states.p_Kv = Y(4); states.m_KCa = Y(5);
    % Ion Concentrations
    states.Ki = Y(6); states.Nai = Y(7); states.Cli = Y(8); states.Caos = Y(9); states.Cais = Y(10); states.Caif = Y(11);
    % Calcium Buffering States
    states.Cab = Y(12); states.Ca_ls = Y(13); states.Ca_hs = Y(14); states.Ca_lf = Y(15); states.Ca_hf = Y(16);
    % Y(17) is a duplicate placeholder for Cab, so it's skipped here.
    % Phototransduction States
    states.cGMP = Y(18); states.PDE_active = Y(19); states.Arr = Y(20); states.Arr_di = Y(21); states.Arr_tetra = Y(22);
    states.G_GTP = Y(23); states.Ga_GDP = Y(24); states.Ga_GTP = Y(25); states.Ga_GTP_PDE_a_Ga_GTP = Y(26);
    states.Ga_GTP_a_PDE_a_Ga_GTP = Y(27); states.Gbg = Y(28); states.Gt = Y(29);
    states.Ops = Y(30); states.Ops_G = Y(31); states.Ops_G_GTP = Y(32); states.Ops_Gt = Y(33);
    states.PDE = Y(34); states.PDE_Ga_GTP = Y(35); states.PDE_a_Ga_GTP = Y(36);
    states.R = Y(37); states.R0 = Y(38); states.R0_G = Y(39); states.R0_G_GTP = Y(40); states.R0_Gt = Y(41);
    states.R0_RKpre = Y(42); states.R1 = Y(43); states.R1_Arr = Y(44); states.R1_G = Y(45); states.R1_G_GTP = Y(46);
    states.R1_Gt = Y(47); states.R1_RKpost = Y(48); states.R1_RKpre = Y(49);
    states.R2 = Y(50); states.R2_Arr = Y(51); states.R2_G = Y(52); states.R2_G_GTP = Y(53); states.R2_Gt = Y(54);
    states.R2_RKpost = Y(55); states.R2_RKpre = Y(56);
    states.R3 = Y(57); states.R3_Arr = Y(58); states.R3_G = Y(59); states.R3_G_GTP = Y(60);
    states.R3_Gt = Y(61); states.R3_RKpost = Y(62); states.R3_RKpre = Y(63);
    states.R4 = Y(64); states.R4_Arr = Y(65); states.R4_G = Y(66); states.R4_G_GTP = Y(67);
    states.R4_Gt = Y(68); states.R4_RKpost = Y(69); states.R4_RKpre = Y(70);
    states.R5 = Y(71); states.R5_Arr = Y(72); states.R5_G = Y(73); states.R5_G_GTP = Y(74);
    states.R5_Gt = Y(75); states.R5_RKpost = Y(76); states.R5_RKpre = Y(77);
    states.R6 = Y(78); states.R6_Arr = Y(79); states.R6_G = Y(80); states.R6_G_GTP = Y(81);
    states.R6_Gt = Y(82); states.R6_RKpost = Y(83); states.R6_RKpre = Y(84);
    states.RGS = Y(85); states.RGS_Ga_GTP_a_PDE_a_Ga_GTP = Y(86); states.RGS_PDE_a_Ga_GTP = Y(87);
    states.RK = Y(88); states.R_Gt = Y(89); states.RecT = Y(90); states.RecR_Ca = Y(91); states.RecR_Ca_RK = Y(92);
    states.Rtot = Y(93); states.E_star = Y(94);
end

function res = unpack_solution_to_results(t, Y, params, Jhv)
    % Converts the ODE solver's output matrix Y into a results struct for plotting.
    n_steps = length(t);
    
    % Pre-allocate memory
    field_names = {'t', 'Vm', 'ICNG', 'INCKX', 'Ih', 'IKv', 'ICaL', 'IKCa', 'IClCa', ...
                   'INaK', 'IPMCA', 'INCX', 'IL_K', 'IL_Na', 'IL_Cl', ...
                   'IL_Caos', 'IL_Cais', 'JNKCC1', 'JKCC2', ...
                   'ATP_Na', 'ATP_Ca', 'ATP_cascade', 'ATP_total', ...
                   'Nai', 'Ki', 'Cli', 'Caos', 'Cais', 'cGMP', 'E_star'};
    for fn = 1:length(field_names)
        res.(field_names{fn}) = zeros(n_steps, 1);
    end
    res.Jhv = Jhv;
    res.t = t;

    % Loop through each time step to calculate and store results
    for k = 1:n_steps
        Y_k = Y(k, :)';
        states_k = unpack_vector_to_states(Y_k);
        currents_k = calculate_currents(states_k, params);
        
        % Calculate cotransporter fluxes (as they are not currents)
        J_NKCC1_k = - params.kNKCC1 * (0.1 * (1 / (1 + exp(16 - params.Ko))) * log(states_k.Ki*states_k.Cli/params.Ko/params.Clo) + log(states_k.Nai*states_k.Cli/params.Nao/params.Clo));
        J_KCC2_k = -params.kKCC2 * (0.3 * log(states_k.Ki*states_k.Cli/params.Ko/params.Clo));

        % Calculate ATP consumption
        stim_k = 0;
        if t(k) >= params.flash_start && t(k) < (params.flash_start + params.flash_duration)
            stim_k = Jhv * 0.43; % 0.43 is the collection area
        end
        [ATP_total_k, ATP_Na_k, ATP_Ca_k, ATP_cascade_k] = calculate_ATP_consumption_accurate(currents_k, states_k, params, stim_k);
        
        % Store state variables
        res.Vm(k) = states_k.Vm;
        res.Ki(k) = states_k.Ki;
        res.Nai(k) = states_k.Nai;
        res.Cli(k) = states_k.Cli;
        res.Caos(k) = states_k.Caos;
        res.Cais(k) = states_k.Cais;
        res.cGMP(k) = states_k.cGMP;
        % Recalculate E_star for storage
        res.E_star(k) = states_k.PDE_a_Ga_GTP + 2 * states_k.Ga_GTP_a_PDE_a_Ga_GTP + states_k.Ga_GTP_PDE_a_Ga_GTP;
        
        % Store currents (in pA)
        res.ICNG(k) = currents_k.I_CNG;
        res.INCKX(k) = currents_k.I_NCKX;
        res.Ih(k) = currents_k.I_h;
        res.IKv(k) = currents_k.I_Kv;
        res.ICaL(k) = currents_k.I_CaL;
        res.IKCa(k) = currents_k.I_KCa;
        res.IClCa(k) = currents_k.I_ClCa;
        res.INaK(k) = currents_k.I_NaK;
        res.INCX(k) = currents_k.I_NCX;
        res.IPMCA(k) = currents_k.I_PMCA;

        % Store leak currents and cotransporter fluxes
        res.IL_K(k) = currents_k.I_L_K;
        res.IL_Na(k) = currents_k.I_L_Na;
        res.IL_Cl(k) = currents_k.I_L_Cl;
        res.IL_Caos(k) = currents_k.I_L_Caos;
        res.IL_Cais(k) = currents_k.I_L_Cais;
        res.JNKCC1(k) = J_NKCC1_k;
        res.JKCC2(k) = J_KCC2_k;

        % Store ATP consumption rates
        res.ATP_Na(k) = ATP_Na_k;
        res.ATP_Ca(k) = ATP_Ca_k;
        res.ATP_cascade(k) = ATP_cascade_k;
        res.ATP_total(k) = ATP_total_k;
    end
end

function analyze_and_print_results(results, t_results, light_intensities, params)
    % Defines time points and variables to output to the console for analysis.
    time_points = [0.5, 0.52, 0.7, 2.0]; % s
    Cm = params.Cm; % pF

    % Define the 24 variables to print, their units, and conversion factors
    vars_to_print = {
        'Vm', 'mV', 1;
        'ICNG', 'pA/pF', 1/Cm; 'INCKX', 'pA/pF', 1/Cm; 'INaK', 'pA/pF', 1/Cm;
        'Ih', 'pA/pF', 1/Cm; 'IKv', 'pA/pF', 1/Cm; 'IKCa', 'pA/pF', 1/Cm;
        'IClCa', 'pA/pF', 1/Cm; 'ICaL', 'pA/pF', 1/Cm; 'IPMCA', 'pA/pF', 1/Cm;
        'INCX', 'fA/pF', 1000/Cm; 'JNKCC1', 'mM/s', 1; 'JKCC2', 'mM/s', 1;
        'IL_K', 'fA/pF', 1000/Cm; 'IL_Na', 'fA/pF', 1000/Cm; 'IL_Cl', 'fA/pF', 1000/Cm;
        'IL_Caos', 'fA/pF', 1000/Cm; 'IL_Cais', 'fA/pF', 1000/Cm;
        'Caos', 'uM', 1; 'Cais', 'uM', 1; 'cGMP', 'uM', 1;
        'Ki', 'mM', 1; 'Nai', 'mM', 1; 'Cli', 'mM', 1;
    };

    for i = 1:length(results)
        fprintf('\n======================================================\n');
        fprintf('  Results Analysis: Light Intensity %.1f photons/um^2/s\n', light_intensities(i));
        fprintf('======================================================\n');
        
        for t_target = time_points
            % Find the index of the data point closest to the target time
            [~, idx] = min(abs(t_results{i} - t_target));
            actual_time = t_results{i}(idx);
            
            fprintf('\n  @ Time Point: %.4f s\n', actual_time);
            fprintf('  ------------------------------------\n');
            
            for v = 1:size(vars_to_print, 1)
                var_name = vars_to_print{v, 1};
                var_unit = vars_to_print{v, 2};
                var_conv = vars_to_print{v, 3};
                
                % Check if the field exists in the results struct
                if isfield(results{i}, var_name)
                    value = results{i}.(var_name)(idx) * var_conv;
                    fprintf('    - %-12s : %12.4f %s\n', var_name, value, var_unit);
                end
            end
        end
    end
end