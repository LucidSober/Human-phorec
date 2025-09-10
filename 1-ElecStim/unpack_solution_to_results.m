%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function res = unpack_solution_to_results(t, Y, Jhv)
    % 将ODE求解器的输出矩阵 Y 转换为绘图函数所需的 results 结构体
    global params states_r

    n_steps = length(t);
    
    % 预分配内存
    field_names = {'t', 'Vm', 'ICNG', 'INCKX', 'Ih', 'IKv', 'ICaL', 'IKCa', 'IClCa', ...
                   'INaK', 'IPMCA', 'INCX', 'IL_K', 'IL_Na', 'IL_Cl', ...
                   'IL_Caos', 'IL_Cais', 'JNKCC1', 'JKCC2', ...
                   'ATP_Na', 'ATP_Ca', 'ATP_cascade', 'ATP_total', ...states_r
                   'Nai', 'Ki', 'Cli', 'Caos', 'Cais', 'cGMP', 'E_star'};
    for fn = 1:length(field_names)
        res.(field_names{fn}) = zeros(n_steps, 1);
    end
    res.Jhv = Jhv;
    res.t = t;

    % 循环遍历每个时间点
    for k = 1:n_steps
        Y_k = Y(k, :);
        states_r = Y2States(Y_k);
        currents_k = cal_ion_currents();
        
        % %%%%%%%%%%%%%%%%%%%%%% 补全代码开始 %%%%%%%%%%%%%%%%%%%%%%%%%%
        % 计算共转运蛋白通量 (因为它们不属于电流)
        J_NKCC1_k = -params.kNKCC1 * (0.1 * (1 / (1 + exp(16 - params.Ko))) * log(states_r.Ki*states_r.Cli/params.Ko/params.Clo) + log(states_r.Nai*states_r.Cli/params.Nao/params.Clo));
        J_KCC2_k = -params.kKCC2 * (0.3 * log(states_r.Ki*states_r.Cli/params.Ko/params.Clo));
        % %%%%%%%%%%%%%%%%%%%%%% 补全代码结束 %%%%%%%%%%%%%%%%%%%%%%%%%%

        % 计算ATP消耗
        stim_k = 0;
        if t(k) >= params.flash_start && t(k) < (params.flash_start + params.flash_duration)
            stim_k = Jhv * 0.43;
        end
        [ATP_Na_k, ATP_Ca_k, ATP_cascade_k] = calculate_ATP_consumption_accurate(currents_k, stim_k);
        
        % 存储状态变量
        res.Vm(k) = states_r.Vm;
        res.Ve(k) = states_r.Ve;
        res.Vi(k) = states_r.Vi;
        res.Ki(k) = states_r.Ki;
        res.Nai(k) = states_r.Nai;
        res.Cli(k) = states_r.Cli; % 之前修正过的
        res.Caos(k) = states_r.Caos;
        res.Cais(k) = states_r.Cais;
        res.cGMP(k) = states_r.cGMP;
        res.E_star(k) = states_r.E_star;
        
        % 存储电流 (单位: pA)
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
        
        % %%%%%%%%%%%%%%%%%%%%%% 补全代码开始 %%%%%%%%%%%%%%%%%%%%%%%%%%
        % 存储漏电流和共转运通量
        res.IL_K(k) = currents_k.I_L_K;
        res.IL_Na(k) = currents_k.I_L_Na;
        res.IL_Cl(k) = currents_k.I_L_Cl;
        res.IL_Caos(k) = currents_k.I_L_Caos;
        res.IL_Cais(k) = currents_k.I_L_Cais;
        res.JNKCC1(k) = J_NKCC1_k;
        res.JKCC2(k) = J_KCC2_k;
        % %%%%%%%%%%%%%%%%%%%%%% 补全代码结束 %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ATP
        res.ATP_Na(k) = ATP_Na_k;
        res.ATP_Ca(k) = ATP_Ca_k;
        res.ATP_cascade(k) = ATP_cascade_k;
        res.ATP_total(k) = ATP_Na_k + ATP_Ca_k + ATP_cascade_k;
    end
end

%[appendix]{"version":"1.0"}
%---
