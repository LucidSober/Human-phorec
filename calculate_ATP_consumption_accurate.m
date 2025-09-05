%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function [ATP_total, ATP_Na, ATP_Ca, ATP_photo] = calculate_ATP_consumption_accurate(currents, stimulus)
% 精确计算ATP消耗
% 输入:
%   currents - 包含所有电流分量的结构体
%   states   - 包含所有状态变量的结构体
%   params   - 包含所有参数的结构体
%   stimulus - 当前的光照强度值
    
global states_r params

    % 提取常量
    F = params.F;
    Vcell = params.V_cell;
    Vis = params.V_is;
    NA = params.NA;
    
    %% 1. ATP for Na+ Pumping
    % 计算总的Na+内流速率 (mM/s)
    Na_influx_currents = currents.I_CNG_Na + 4 * currents.I_NCKX + currents.I_h_Na + ...
                         currents.I_CaL_Na + 3 * currents.I_NCX + currents.I_L_Na;
    
    Na_influx_rate = -Na_influx_currents * 1e-12 / (F * Vcell) * 1e3; % mol/s -> mM/s
    
    % INaK泵出3个Na+消耗1个ATP (molecules/s)
    ATP_Na = abs(Na_influx_rate) * 1e-3 * NA * Vcell / 3;
    
    %% 2. ATP for Ca2+ Pumping
    % 计算内段的总Ca2+内流速率 (uM/s)
    Ca_influx_currents = currents.I_CaL_Ca + currents.I_L_Cais - 2 * currents.I_NCX;
    
    Ca_influx_rate = -Ca_influx_currents * 1e-12 / (2 * F * Vis) * 1e6; % mol/s -> uM/s
    
    % IPMCA泵出1个Ca2+消耗1个ATP (molecules/s)
    ATP_Ca = abs(Ca_influx_rate) * 1e-6 * NA * Vis;
    
    %% 3. 光转导级联的ATP消耗
    ATP_photo = 0; % 默认为0
    if stimulus > 0
        % 在光下，ATP主要消耗在视紫红质磷酸化和G蛋白激活上
        % 为此，我们需要重新计算与ATP消耗相关的反应速率
        v_r3_0 = params.kRK3_ATP * states_r.R0_RKpre;
        v_r3_1 = params.kRK3_ATP * states_r.R1_RKpre;
        v_r3_2 = params.kRK3_ATP * states_r.R2_RKpre;
        v_r3_3 = params.kRK3_ATP * states_r.R3_RKpre;
        v_r3_4 = params.kRK3_ATP * states_r.R4_RKpre;
        v_r3_5 = params.kRK3_ATP * states_r.R5_RKpre;
        v_r10 = params.kG5_GTP * states_r.Ops_G;
        v_r15_0 = params.kG5_GTP * states_r.R0_G;
        v_r15_1 = params.kG5_GTP * states_r.R1_G;
        v_r15_2 = params.kG5_GTP * states_r.R2_G;
        v_r15_3 = params.kG5_GTP * states_r.R3_G;
        v_r15_4 = params.kG5_GTP * states_r.R4_G;
        v_r15_5 = params.kG5_GTP * states_r.R5_G;
        v_r15_6 = params.kG5_GTP * states_r.R6_G;
        
        % 将反应速率(uM/s)转换为分子数/秒
        ATP_photo_rate_uMs = v_r3_0 + v_r3_1 + v_r3_2 + v_r3_3 + v_r3_4 + v_r3_5 + ...
                             v_r10 + v_r15_0 + v_r15_1 + v_r15_2 + v_r15_3 + v_r15_4 + ...
                             v_r15_5 + v_r15_6;
        ATP_photo = ATP_photo_rate_uMs * 1e-6 * params.Vcyto * 1e-12 * NA;
    end
    
    %% 4. 总ATP消耗
    ATP_total = ATP_Na + ATP_Ca + ATP_photo;
end

%[appendix]{"version":"1.0"}
%---
