% main.m
% This script runs the complete rod photoreceptor model simulation with a single, specific stimulus.

clear all; close all; clc;

%% 1. Initialization
fprintf('1. Initializing parameters and state variables...\n');
params = init_parameters();      
states = init_state_variables();  
Y0 = pack_states_to_vector(states); 

%% 2. Simulation Setup
% ODE solver options - using stricter tolerances for accuracy
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, ...
                 'MaxStep', 0.001, ...
                 'InitialStep', 1e-6);

%% 3. Find the Steady State in Darkness
fprintf('2. Finding the steady state in dark conditions...\n');
tspan_ss = [0 10]; % Increased duration for stability
[~, Y_ss] = ode15s(@(t,Y) derivatives(t, Y, params), tspan_ss, Y0, options);
Y_steady_state = Y_ss(end, :)';

states_steady_state = unpack_vector_to_states(Y_steady_state);
fprintf('   Steady state achieved!\n\n');

%% 4. Run SINGLE SIMULATION with the new Gaussian Stimulus
fprintf('3. Starting simulation with the new, low-amplitude Gaussian stimulus...\n');
t_total = 5.0; % Total simulation time in seconds
[t_full, Y_full] = ode15s(@(t,Y) derivatives(t, Y, params), [0, t_total], Y_steady_state, options);

% Unpack the full solution into a more accessible format
results = unpack_solution_to_results(t_full, Y_full, params);
fprintf('   Simulation complete.\n\n');

%% 5. Plotting Results
fprintf('4. Plotting results...\n');

% --- Plot 1: Stimulus Profile ---
figure('Name', 'Stimulus Profile');
plot(results.t, results.stimulus, 'r-', 'LineWidth', 2);
title('Applied Light Stimulus (Low-Amplitude Gaussian Waveform)');
xlabel('Time (s)');
ylabel('Scaled Stimulus Intensity');
grid on;
box on;
ax = gca;
ax.YAxis.Exponent = 0; % Prevent scientific notation on y-axis

% --- Plot 2: Key Model Responses ---
figure('Name', 'Model Responses');

% Membrane Potential
subplot(2, 2, 1);
plot(results.t, results.Vm, 'k-', 'LineWidth', 1.5);
title('Membrane Potential (Vm)');
xlabel('Time (s)');
ylabel('Voltage (mV)');
grid on;

% cGMP Concentration
subplot(2, 2, 2);
plot(results.t, results.cGMP, 'b-', 'LineWidth', 1.5);
title('cGMP Concentration');
xlabel('Time (s)');
ylabel('Concentration (μM)');
grid on;

% Outer Segment Calcium
subplot(2, 2, 3);
plot(results.t, results.Caos, 'g-', 'LineWidth', 1.5);
title('Outer Segment Ca2+ (Caos)');
xlabel('Time (s)');
ylabel('Concentration (μM)');
grid on;

% Total Circulating Current (Dark Current)
subplot(2, 2, 4);
total_current = results.ICNG + results.INCKX;
plot(results.t, total_current, 'm-', 'LineWidth', 1.5);
title('Total Circulating Current (ICNG + INCKX)');
xlabel('Time (s)');
ylabel('Current (pA)');
grid on;

fprintf('Simulations and plotting complete!\n');

%% Helper Functions

function Y0 = pack_states_to_vector(states)
    Y0 = zeros(94, 1);
    Y0(1) = states.Vm; Y0(2) = states.p_CaL; Y0(3) = states.p_h; Y0(4) = states.p_Kv; Y0(5) = states.m_KCa;
    Y0(6) = states.Ki; Y0(7) = states.Nai; Y0(8) = states.Cli; Y0(9) = states.Caos; Y0(10) = states.Cais; Y0(11) = states.Caif;
    Y0(12) = states.Cab; Y0(13) = states.Ca_ls; Y0(14) = states.Ca_hs; Y0(15) = states.Ca_lf; Y0(16) = states.Ca_hf; Y0(17) = states.Cab;
    Y0(18) = states.cGMP; Y0(19) = states.PDE_active;
    Y0(20) = states.Arr; Y0(21) = states.Arr_di; Y0(22) = states.Arr_tetra; Y0(23) = states.G_GTP; Y0(24) = states.Ga_GDP; Y0(25) = states.Ga_GTP; Y0(26) = states.Ga_GTP_PDE_a_Ga_GTP; Y0(27) = states.Ga_GTP_a_PDE_a_Ga_GTP; Y0(28) = states.Gbg; Y0(29) = states.Gt; Y0(30) = states.Ops; Y0(31) = states.Ops_G; Y0(32) = states.Ops_G_GTP; Y0(33) = states.Ops_Gt; Y0(34) = states.PDE; Y0(35) = states.PDE_Ga_GTP; Y0(36) = states.PDE_a_Ga_GTP; Y0(37) = states.R; Y0(38) = states.R0; Y0(39) = states.R0_G; Y0(40) = states.R0_G_GTP; Y0(41) = states.R0_Gt; Y0(42) = states.R0_RKpre; Y0(43) = states.R1; Y0(44) = states.R1_Arr; Y0(45) = states.R1_G; Y0(46) = states.R1_G_GTP; Y0(47) = states.R1_Gt; Y0(48) = states.R1_RKpost; Y0(49) = states.R1_RKpre; Y0(50) = states.R2; Y0(51) = states.R2_Arr; Y0(52) = states.R2_G; Y0(53) = states.R2_G_GTP; Y0(54) = states.R2_Gt; Y0(55) = states.R2_RKpost; Y0(56) = states.R2_RKpre; Y0(57) = states.R3; Y0(58) = states.R3_Arr; Y0(59) = states.R3_G; Y0(60) = states.R3_G_GTP; Y0(61) = states.R3_Gt; Y0(62) = states.R3_RKpost; Y0(63) = states.R3_RKpre; Y0(64) = states.R4; Y0(65) = states.R4_Arr; Y0(66) = states.R4_G; Y0(67) = states.R4_G_GTP; Y0(68) = states.R4_Gt; Y0(69) = states.R4_RKpost; Y0(70) = states.R4_RKpre; Y0(71) = states.R5; Y0(72) = states.R5_Arr; Y0(73) = states.R5_G; Y0(74) = states.R5_G_GTP; Y0(75) = states.R5_Gt; Y0(76) = states.R5_RKpost; Y0(77) = states.R5_RKpre; Y0(78) = states.R6; Y0(79) = states.R6_Arr; Y0(80) = states.R6_G; Y0(81) = states.R6_G_GTP; Y0(82) = states.R6_Gt; Y0(83) = states.R6_RKpost; Y0(84) = states.R6_RKpre; Y0(85) = states.RGS; Y0(86) = states.RGS_Ga_GTP_a_PDE_a_Ga_GTP; Y0(87) = states.RGS_PDE_a_Ga_GTP; Y0(88) = states.RK; Y0(89) = states.R_Gt; Y0(90) = states.RecT; Y0(91) = states.RecR_Ca; Y0(92) = states.RecR_Ca_RK; Y0(93) = states.Rtot; Y0(94) = states.E_star;
end

function states = unpack_vector_to_states(Y)
    states = struct();
    states.Vm = Y(1); states.p_CaL = Y(2); states.p_h = Y(3); states.p_Kv = Y(4); states.m_KCa = Y(5);
    states.Ki = Y(6); states.Nai = Y(7); states.Cli = Y(8); states.Caos = Y(9); states.Cais = Y(10); states.Caif = Y(11);
    states.Cab = Y(12); states.Ca_ls = Y(13); states.Ca_hs = Y(14); states.Ca_lf = Y(15); states.Ca_hf = Y(16);
    states.cGMP = Y(18); states.PDE_active = Y(19); states.Arr = Y(20); states.Arr_di = Y(21); states.Arr_tetra = Y(22); states.G_GTP = Y(23); states.Ga_GDP = Y(24); states.Ga_GTP = Y(25); states.Ga_GTP_PDE_a_Ga_GTP = Y(26); states.Ga_GTP_a_PDE_a_Ga_GTP = Y(27); states.Gbg = Y(28); states.Gt = Y(29); states.Ops = Y(30); states.Ops_G = Y(31); states.Ops_G_GTP = Y(32); states.Ops_Gt = Y(33); states.PDE = Y(34); states.PDE_Ga_GTP = Y(35); states.PDE_a_Ga_GTP = Y(36); states.R = Y(37); states.R0 = Y(38); states.R0_G = Y(39); states.R0_G_GTP = Y(40); states.R0_Gt = Y(41); states.R0_RKpre = Y(42); states.R1 = Y(43); states.R1_Arr = Y(44); states.R1_G = Y(45); states.R1_G_GTP = Y(46); states.R1_Gt = Y(47); states.R1_RKpost = Y(48); states.R1_RKpre = Y(49); states.R2 = Y(50); states.R2_Arr = Y(51); states.R2_G = Y(52); states.R2_G_GTP = Y(53); states.R2_Gt = Y(54); states.R2_RKpost = Y(55); states.R2_RKpre = Y(56); states.R3 = Y(57); states.R3_Arr = Y(58); states.R3_G = Y(59); states.R3_G_GTP = Y(60); states.R3_Gt = Y(61); states.R3_RKpost = Y(62); states.R3_RKpre = Y(63); states.R4 = Y(64); states.R4_Arr = Y(65); states.R4_G = Y(66); states.R4_G_GTP = Y(67); states.R4_Gt = Y(68); states.R4_RKpost = Y(69); states.R4_RKpre = Y(70); states.R5 = Y(71); states.R5_Arr = Y(72); states.R5_G = Y(73); states.R5_G_GTP = Y(74); states.R5_Gt = Y(75); states.R5_RKpost = Y(76); states.R5_RKpre = Y(77); states.R6 = Y(78); states.R6_Arr = Y(79); states.R6_G = Y(80); states.R6_G_GTP = Y(81); states.R6_Gt = Y(82); states.R6_RKpost = Y(83); states.R6_RKpre = Y(84); states.RGS = Y(85); states.RGS_Ga_GTP_a_PDE_a_Ga_GTP = Y(86); states.RGS_PDE_a_Ga_GTP = Y(87); states.RK = Y(88); states.R_Gt = Y(89); states.RecT = Y(90); states.RecR_Ca = Y(91); states.RecR_Ca_RK = Y(92); states.Rtot = Y(93); states.E_star = Y(94);
end

function res = unpack_solution_to_results(t, Y, params)
    n_steps = length(t);
    field_names = {'t', 'Vm', 'ICNG', 'INCKX', 'Ih', 'IKv', 'ICaL', 'IKCa', 'IClCa', ...
                   'INaK', 'IPMCA', 'INCX', 'IL_K', 'IL_Na', 'IL_Cl', ...
                   'IL_Caos', 'IL_Cais', 'JNKCC1', 'JKCC2', ...
                   'ATP_Na', 'ATP_Ca', 'ATP_cascade', 'ATP_total', ...
                   'Nai', 'Ki', 'Cli', 'Caos', 'Cais', 'cGMP', 'E_star', 'stimulus'};
    for fn = 1:length(field_names)
        res.(field_names{fn}) = zeros(n_steps, 1);
    end
    res.t = t;

    % Define Gaussian stimulus parameters here for consistency
    t_peak = 1.75;
    tau_decay_smoother = 0.6;
    sigma = tau_decay_smoother / sqrt(2);
    % Use the same, new peak intensity as in derivatives.m
    peak_Jhv = 210;
    max_stim_model = peak_Jhv * 0.43;

    for k = 1:n_steps
        Y_k = Y(k, :)';
        states_k = unpack_vector_to_states(Y_k);
        currents_k = calculate_currents(states_k, params);
        
        % Calculate and store the stimulus profile at time t(k)
        t_k = t(k);
        stim_waveform = exp( - (t_k - t_peak)^2 / (2 * sigma^2) );
        stim_k = stim_waveform * max_stim_model;
        res.stimulus(k) = stim_k;
        
        % Use this consistent stimulus value for ATP calculation
        [~, ATP_Na_k, ATP_Ca_k, ATP_cascade_k] = calculate_ATP_consumption_accurate(currents_k, states_k, params, stim_k);
        
        % Store state variables and currents
        res.Vm(k) = states_k.Vm;
        res.Ki(k) = states_k.Ki;
        res.Nai(k) = states_k.Nai;
        res.Cli(k) = states_k.Cli;
        res.Caos(k) = states_k.Caos;
        res.Cais(k) = states_k.Cais;
        res.cGMP(k) = states_k.cGMP;
        res.ICNG(k) = currents_k.I_CNG;
        res.INCKX(k) = currents_k.I_NCKX;
        
        % Store ATP consumption rates
        res.ATP_Na(k) = ATP_Na_k;
        res.ATP_Ca(k) = ATP_Ca_k;
        res.ATP_cascade(k) = ATP_cascade_k;
        res.ATP_total(k) = ATP_Na_k + ATP_Ca_k + ATP_cascade_k;
    end
end