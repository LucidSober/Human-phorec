% main.m (Final Version - No Global, with Ve(t) check)

clear all; close all; clc;

%% ========================================================================
%  PART 1: 运行健康模型，生成并保存目标响应曲线
% =========================================================================
fprintf('PART 1: Generating target response from a healthy cell...\n');

% 1.1 初始化
params_healthy = init_parameters();
params_healthy.flash_start = 0;
params_healthy.flash_duration = 0;
params_healthy.Jhv = 0;
states_healthy = init_state_variables();
Y0_healthy = pack_states_to_vector(states_healthy);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.001);

% 1.2 为健康模型定义 Ve_func (始终为0)
Ve_func_healthy = @(t) 0 * t;

% 1.3 寻找暗适应稳态
fprintf('  Finding steady state in darkness...\n');
[~, Y_ss_h] = ode15s(@(t,Y) derivatives(t, Y, params_healthy, Ve_func_healthy), [0 5], Y0_healthy, options);
Y_steady_state_h = Y_ss_h(end, :)';

% 1.4 运行闪光模拟
fprintf('  Simulating light flash response...\n');
light_intensities = [1.7, 4.8, 15.2, 39.4, 125, 444, 1406, 4630];
params_healthy.Jhv = light_intensities(4);
params_healthy.flash_start = 0.5;
params_healthy.flash_duration = 0.02;
[t_healthy, Y_healthy] = ode15s(@(t,Y) derivatives(t, Y, params_healthy, Ve_func_healthy), [0, 5.0], Y_steady_state_h, options);
Vm_healthy = Y_healthy(:,1);

% 1.5 保存目标响应
fprintf('  Saving healthy response data...\n\n');
save('target_response.mat', 't_healthy', 'Vm_healthy');


%% ========================================================================
%  PART 2: 运行病变模型，并施加外部电刺激
% =========================================================================
fprintf('PART 2: Simulating the "diseased" cell with electrical stimulation...\n');

% 2.1 初始化
params_diseased = init_parameters();
params_diseased.flash_start = 0;
params_diseased.flash_duration = 0;
params_diseased.Jhv = 0;
states_diseased = init_state_variables();
Y0_diseased = pack_states_to_vector(states_diseased);
Ve_func_zero = @(t) 0 * t;
[~, Y_ss_d] = ode15s(@(t,Y) derivatives(t, Y, params_diseased, Ve_func_zero), [0 5], Y0_diseased, options);
Y_steady_state_d = Y_ss_d(end, :)';

% 2.2 设置电刺激参数
fprintf('  Setting up extracellular stimulation...\n');
params_diseased.Rf = 200; params_diseased.Cdl = 1e-6; params_diseased.Rs = 100;
Tp = 0.5; I_amp = 35e-3; delay = 0.5;
t_total_stim = 5.0; dt = 1e-5;
t_stim = 0:dt:t_total_stim;
I_stim = zeros(size(t_stim));
I_stim(t_stim >= delay & t_stim < (delay + Tp)) = I_amp;

% 2.3 计算细胞外电压 Ve(t)
fprintf('  Calculating extracellular voltage Ve(t)...\n');
I_func = @(t) interp1(t_stim, I_stim, t, 'previous', 0);
dVc_dt = @(t, Vc) (I_func(t) - Vc/params_diseased.Rf) / params_diseased.Cdl;
[t_out_Ve, Vc_out] = ode15s(dVc_dt, t_stim, 0);
Ve_out = Vc_out + I_func(t_out_Ve) * params_diseased.Rs;
Ve_func_diseased = @(te) interp1(t_out_Ve, Ve_out, te, 'linear', 'extrap');

% 2.4 [调试步骤] 检查生成的 Ve(t) 信号
figure('Name', 'DEBUG: Check of Generated Ve(t) Signal');
plot(t_out_Ve * 1000, Ve_out, 'k-', 'LineWidth', 2);
title('Generated V_e(t) Signal');
xlabel('Time (ms)'); ylabel('V_e (mV)'); grid on;
fprintf('  >>> A plot window named "DEBUG" has been generated.\n');
fprintf('  >>> Please CHECK if the Ve(t) signal in this plot is correct (not a flat line).\n');
fprintf('  >>> Press any key to continue the main simulation...\n');
%pause;

% 2.5 运行主仿真
fprintf('  Running main ODE solver for the stimulated cell...\n');
[t_diseased, Y_diseased] = ode15s(@(t,Y) derivatives(t, Y, params_diseased, Ve_func_diseased), [0, t_total_stim], Y_steady_state_d, options);
Vi_diseased = Y_diseased(:,1);
Vm_stimulated = Vi_diseased - Ve_func_diseased(t_diseased);
fprintf('Simulations complete.\n\n');

%% ========================================================================
%  PART 3: 结果可视化与比较
% =========================================================================
fprintf('PART 3: Plotting final results...\n');
load('target_response.mat');
figure('Name', 'Photoreceptor Stimulation Comparison', 'Position', [100, 100, 800, 600]);
hold on;
plot(t_healthy * 1000, Vm_healthy, 'g--', 'LineWidth', 3, 'DisplayName', 'Healthy Response to Light (Target)');
plot(t_diseased * 1000, Vm_stimulated, 'b-', 'LineWidth', 2, 'DisplayName', 'Stimulated "Diseased" Cell Response');
hold off;
title('Comparison of Membrane Potential V_m(t)');
xlabel('Time (ms)'); ylabel('V_m (mV)');
legend('show', 'Location', 'southeast'); grid on; box on;
ylim([-65, -30]);


% (以下是原始代码中的辅助函数，需要保留在文件末尾或作为单独文件)
function Y0 = pack_states_to_vector(states)
    Y0 = zeros(94, 1);
    Y0(1) = states.Vm; Y0(2) = states.p_CaL; Y0(3) = states.p_h; Y0(4) = states.p_Kv; Y0(5) = states.m_KCa;
    Y0(6) = states.Ki; Y0(7) = states.Nai; Y0(8) = states.Cli; Y0(9) = states.Caos; Y0(10) = states.Cais; Y0(11) = states.Caif;
    Y0(12) = states.Cab; Y0(13) = states.Ca_ls; Y0(14) = states.Ca_hs; Y0(15) = states.Ca_lf; Y0(16) = states.Ca_hf; Y0(17) = states.Cab;
    Y0(18) = states.cGMP; Y0(19) = states.PDE_active;
    Y0(20) = states.Arr; Y0(21) = states.Arr_di; Y0(22) = states.Arr_tetra; Y0(23) = states.G_GTP; Y0(24) = states.Ga_GDP; Y0(25) = states.Ga_GTP;
    Y0(26) = states.Ga_GTP_PDE_a_Ga_GTP; Y0(27) = states.Ga_GTP_a_PDE_a_Ga_GTP; Y0(28) = states.Gbg; Y0(29) = states.Gt;
    Y0(30) = states.Ops; Y0(31) = states.Ops_G; Y0(32) = states.Ops_G_GTP; Y0(33) = states.Ops_Gt;
    Y0(34) = states.PDE; Y0(35) = states.PDE_Ga_GTP; Y0(36) = states.PDE_a_Ga_GTP; Y0(37) = states.R; Y0(38) = states.R0; Y0(39) = states.R0_G; Y0(40) = states.R0_G_GTP;
    Y0(41) = states.R0_Gt; Y0(42) = states.R0_RKpre; Y0(43) = states.R1; Y0(44) = states.R1_Arr; Y0(45) = states.R1_G; Y0(46) = states.R1_G_GTP;
    Y0(47) = states.R1_Gt; Y0(48) = states.R1_RKpost; Y0(49) = states.R1_RKpre; Y0(50) = states.R2; Y0(51) = states.R2_Arr; Y0(52) = states.R2_G; Y0(53) = states.R2_G_GTP;
    Y0(54) = states.R2_Gt; Y0(55) = states.R2_RKpost; Y0(56) = states.R2_RKpre; Y0(57) = states.R3; Y0(58) = states.R3_Arr; Y0(59) = states.R3_G; Y0(60) = states.R3_G_GTP;
    Y0(61) = states.R3_Gt; Y0(62) = states.R3_RKpost; Y0(63) = states.R3_RKpre; Y0(64) = states.R4; Y0(65) = states.R4_Arr; Y0(66) = states.R4_G; Y0(67) = states.R4_G_GTP;
    Y0(68) = states.R4_Gt; Y0(69) = states.R4_RKpost; Y0(70) = states.R4_RKpre; Y0(71) = states.R5; Y0(72) = states.R5_Arr; Y0(73) = states.R5_G; Y0(74) = states.R5_G_GTP;
    Y0(75) = states.R5_Gt; Y0(76) = states.R5_RKpost; Y0(77) = states.R5_RKpre; Y0(78) = states.R6; Y0(79) = states.R6_Arr; Y0(80) = states.R6_G; Y0(81) = states.R6_G_GTP;
    Y0(82) = states.R6_Gt; Y0(83) = states.R6_RKpost; Y0(84) = states.R6_RKpre; Y0(85) = states.RGS; Y0(86) = states.RGS_Ga_GTP_a_PDE_a_Ga_GTP; Y0(87) = states.RGS_PDE_a_Ga_GTP;
    Y0(88) = states.RK; Y0(89) = states.R_Gt; Y0(90) = states.RecT; Y0(91) = states.RecR_Ca; Y0(92) = states.RecR_Ca_RK;
    Y0(93) = states.Rtot; Y0(94) = states.E_star;
end
