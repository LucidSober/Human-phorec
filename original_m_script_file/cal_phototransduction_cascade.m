function [dStates, E_star, v] = cal_phototransduction_cascade()
% 计算光转导级联反应中所有状态变量的导数 (d/dt)
%
% 输出:
%   dStates   - 包含所有状态变量导数(d/dt)的结构体
%   E_star    - 活化PDE的浓度
%   v         - 包含所有反应速率的结构体
  
global params states_r

%% 1. 计算速率常数
    kRK1_1 = params.kRK1_0 * exp(-params.omega);
    kRK1_2 = params.kRK1_0 * exp(-params.omega * 2);
    kRK1_3 = params.kRK1_0 * exp(-params.omega * 3);
    kRK1_4 = params.kRK1_0 * exp(-params.omega * 4);
    kRK1_5 = params.kRK1_0 * exp(-params.omega * 5);
    kRK1_6 = 0;
    
    kA1_1 = params.kArr;
    kA1_2 = params.kArr + 1 * params.m_Arr;
    kA1_3 = params.kArr + 2 * params.m_Arr;
    kA1_4 = params.kArr + 3 * params.m_Arr;
    kA1_5 = params.kArr + 3 * params.m_Arr;
    kA1_6 = params.kArr + 3 * params.m_Arr;
    
    kGpre1 = params.kG1_0 * 1.6;
    kGpre2 = params.kG2 * 315;
    kG1_1 = params.kG1_0 * exp(-params.omega_G);
    kG1_2 = params.kG1_0 * exp(-params.omega_G * 2);
    kG1_3 = params.kG1_0 * exp(-params.omega_G * 3);
    kG1_4 = params.kG1_0 * exp(-params.omega_G * 4);
    kG1_5 = params.kG1_0 * exp(-params.omega_G * 5);
    kG1_6 = params.kG1_0 * exp(-params.omega_G * 6);
    
    %% 2. 计算反应速率 (v)
    v.r1 = params.stimulus * states_r.R / params.Rtot;
    v.rstprec = params.stimulus * states_r.R_Gt / params.Rtot;
    
    % R-RK 反应
    v.r2_0 = params.kRK1_0 * states_r.RK * states_r.R0 - params.kRK2 * states_r.R0_RKpre;
    v.r2_1 = kRK1_1 * states_r.RK * states_r.R1 - params.kRK2 * states_r.R1_RKpre;
    v.r2_2 = kRK1_2 * states_r.RK * states_r.R2 - params.kRK2 * states_r.R2_RKpre;
    v.r2_3 = kRK1_3 * states_r.RK * states_r.R3 - params.kRK2 * states_r.R3_RKpre;
    v.r2_4 = kRK1_4 * states_r.RK * states_r.R4 - params.kRK2 * states_r.R4_RKpre;
    v.r2_5 = kRK1_5 * states_r.RK * states_r.R5 - params.kRK2 * states_r.R5_RKpre;
    v.r2_6 = kRK1_6 * states_r.RK * states_r.R6 - params.kRK2 * states_r.R6_RKpre;
    
    v.r3_0 = params.kRK3_ATP * states_r.R0_RKpre;
    v.r3_1 = params.kRK3_ATP * states_r.R1_RKpre;
    v.r3_2 = params.kRK3_ATP * states_r.R2_RKpre;
    v.r3_3 = params.kRK3_ATP * states_r.R3_RKpre;
    v.r3_4 = params.kRK3_ATP * states_r.R4_RKpre;
    v.r3_5 = params.kRK3_ATP * states_r.R5_RKpre;
    
    v.r4_1 = params.kRK4 * states_r.R1_RKpost;
    v.r4_2 = params.kRK4 * states_r.R2_RKpost;
    v.r4_3 = params.kRK4 * states_r.R3_RKpost;
    v.r4_4 = params.kRK4 * states_r.R4_RKpost;
    v.r4_5 = params.kRK4 * states_r.R5_RKpost;
    v.r4_6 = params.kRK4 * states_r.R6_RKpost;
    
    % Arrestin 反应
    v.r5_1 = kA1_1 * states_r.Arr * states_r.R1 - params.kA2 * states_r.R1_Arr;
    v.r5_2 = kA1_2 * states_r.Arr * states_r.R2 - params.kA2 * states_r.R2_Arr;
    v.r5_3 = kA1_3 * states_r.Arr * states_r.R3 - params.kA2 * states_r.R3_Arr;
    v.r5_4 = kA1_4 * states_r.Arr * states_r.R4 - params.kA2 * states_r.R4_Arr;
    v.r5_5 = kA1_5 * states_r.Arr * states_r.R5 - params.kA2 * states_r.R5_Arr;
    v.r5_6 = kA1_6 * states_r.Arr * states_r.R6 - params.kA2 * states_r.R6_Arr;
    
    v.r6_1 = params.kA3 * states_r.R1_Arr;
    v.r6_2 = params.kA3 * states_r.R2_Arr;
    v.r6_3 = params.kA3 * states_r.R3_Arr;
    v.r6_4 = params.kA3 * states_r.R4_Arr;
    v.r6_5 = params.kA3 * states_r.R5_Arr;
    v.r6_6 = params.kA3 * states_r.R6_Arr;
    
    % 热激活
    v.r7_0 = params.ktherm * states_r.R0;
    v.r7_1 = params.ktherm * states_r.R1;
    v.r7_2 = params.ktherm * states_r.R2;
    v.r7_3 = params.ktherm * states_r.R3;
    v.r7_4 = params.ktherm * states_r.R4;
    v.r7_5 = params.ktherm * states_r.R5;
    v.r7_6 = params.ktherm * states_r.R6;
    
    % G蛋白反应
    v.r8 = params.kOps * states_r.Ops * states_r.Gt - params.kG2 * states_r.Ops_Gt;
    v.r9 = params.kG3 * states_r.Ops_Gt - params.kG4_GDP * states_r.Ops_G;
    v.r10 = params.kG5_GTP * states_r.Ops_G;
    v.r11 = params.kG6 * states_r.Ops_G_GTP;
    v.r12 = params.kRrecyc * states_r.Ops;
    
    v.GtRpre = kGpre1 * states_r.Gt * states_r.R - kGpre2 * states_r.R_Gt;
    
    v.r13_0 = params.kG1_0 * states_r.Gt * states_r.R0 - params.kG2 * states_r.R0_Gt;
    v.r13_1 = kG1_1 * states_r.Gt * states_r.R1 - params.kG2 * states_r.R1_Gt;
    v.r13_2 = kG1_2 * states_r.Gt * states_r.R2 - params.kG2 * states_r.R2_Gt;
    v.r13_3 = kG1_3 * states_r.Gt * states_r.R3 - params.kG2 * states_r.R3_Gt;
    v.r13_4 = kG1_4 * states_r.Gt * states_r.R4 - params.kG2 * states_r.R4_Gt;
    v.r13_5 = kG1_5 * states_r.Gt * states_r.R5 - params.kG2 * states_r.R5_Gt;
    v.r13_6 = kG1_6 * states_r.Gt * states_r.R6 - params.kG2 * states_r.R6_Gt;
    
    v.r14_0 = params.kG3 * states_r.R0_Gt - params.kG4_GDP * states_r.R0_G;
    v.r14_1 = params.kG3 * states_r.R1_Gt - params.kG4_GDP * states_r.R1_G;
    v.r14_2 = params.kG3 * states_r.R2_Gt - params.kG4_GDP * states_r.R2_G;
    v.r14_3 = params.kG3 * states_r.R3_Gt - params.kG4_GDP * states_r.R3_G;
    v.r14_4 = params.kG3 * states_r.R4_Gt - params.kG4_GDP * states_r.R4_G;
    v.r14_5 = params.kG3 * states_r.R5_Gt - params.kG4_GDP * states_r.R5_G;
    v.r14_6 = params.kG3 * states_r.R6_Gt - params.kG4_GDP * states_r.R6_G;
    
    v.r15_0 = params.kG5_GTP * states_r.R0_G;
    v.r15_1 = params.kG5_GTP * states_r.R1_G;
    v.r15_2 = params.kG5_GTP * states_r.R2_G;
    v.r15_3 = params.kG5_GTP * states_r.R3_G;
    v.r15_4 = params.kG5_GTP * states_r.R4_G;
    v.r15_5 = params.kG5_GTP * states_r.R5_G;
    v.r15_6 = params.kG5_GTP * states_r.R6_G;
    
    v.r16_0 = params.kG6 * states_r.R0_G_GTP;
    v.r16_1 = params.kG6 * states_r.R1_G_GTP;
    v.r16_2 = params.kG6 * states_r.R2_G_GTP;
    v.r16_3 = params.kG6 * states_r.R3_G_GTP;
    v.r16_4 = params.kG6 * states_r.R4_G_GTP;
    v.r16_5 = params.kG6 * states_r.R5_G_GTP;
    v.r16_6 = params.kG6 * states_r.R6_G_GTP;
    
    v.r17 = params.kG7 * states_r.G_GTP;
    
    % PDE反应
    v.r18 = params.kP1 * states_r.PDE * states_r.Ga_GTP - params.kP1_rev * states_r.PDE_Ga_GTP;
    v.r19 = params.kP2 * states_r.PDE_Ga_GTP;
    v.r20 = params.kP3 * states_r.PDE_a_Ga_GTP * states_r.Ga_GTP;
    v.r21 = params.kP4 * states_r.Ga_GTP_PDE_a_Ga_GTP;
    
    % RGS反应
    v.r22 = params.kRGS1 * states_r.RGS * states_r.Ga_GTP_a_PDE_a_Ga_GTP;
    v.r23 = params.kRGS2 * states_r.RGS_Ga_GTP_a_PDE_a_Ga_GTP;
    v.r24 = params.kRGS1 * states_r.RGS * states_r.PDE_a_Ga_GTP;
    v.r25 = params.kRGS2 * states_r.RGS_PDE_a_Ga_GTP;
    
    v.r26 = params.kPDEshutoff * states_r.PDE_a_Ga_GTP;
    v.r27 = params.kPDEshutoff * states_r.Ga_GTP_a_PDE_a_Ga_GTP;
    v.r28 = params.kGshutoff * states_r.Ga_GTP;
    v.r29 = params.kGrecyc * states_r.Gbg * states_r.Ga_GDP;
    
    % Recoverin反应
    v.r30 = params.kRec1 * states_r.RecT * states_r.Caos - params.kRec2 * states_r.RecR_Ca;
    v.r31 = params.kRec3 * states_r.RecR_Ca * states_r.RK - params.kRec4 * states_r.RecR_Ca_RK;
    
    % Arrestin二聚化
    v.r_diarr = params.kA4 * states_r.Arr * states_r.Arr - params.kA5 * states_r.Arr_di;
    v.r_tetraarr = params.kA4 * states_r.Arr_di * states_r.Arr_di - params.kA5 * states_r.Arr_tetra;
    
    %% 3. 计算E_star
    E_star = states_r.PDE_a_Ga_GTP + 2 * states_r.Ga_GTP_a_PDE_a_Ga_GTP + states_r.Ga_GTP_PDE_a_Ga_GTP;
    
    %% 4. 计算状态导数 (d/dt)
    dStates = struct();
    
    % 所有状态变量的导数... (根据 Invergo.txt)
    dStates.Arr = (-v.r5_1 - v.r5_2 - v.r5_3 - v.r5_4 - v.r5_5 - v.r5_6 + v.r6_1 + v.r6_2 + v.r6_3 + v.r6_4 + v.r6_5 + v.r6_6 - v.r_diarr - v.r_diarr);
    dStates.Arr_di = v.r_diarr - v.r_tetraarr - v.r_tetraarr;
    dStates.Arr_tetra = v.r_tetraarr;
    dStates.G_GTP = v.r11 + v.r16_0 + v.r16_1 + v.r16_2 + v.r16_3 + v.r16_4 + v.r16_5 + v.r16_6 - v.r17;
    dStates.Ga_GDP = v.r23 + v.r25 + v.r26 + v.r27 + v.r28 - v.r29;
    dStates.Ga_GTP = v.r17 - v.r18 - v.r20 - v.r28;
    dStates.Ga_GTP_PDE_a_Ga_GTP = v.r20 - v.r21;
    dStates.Ga_GTP_a_PDE_a_Ga_GTP = v.r21 - v.r22 - v.r27;
    dStates.Gbg = v.r17 - v.r29;
    dStates.Gt = -v.r8 - v.GtRpre - v.r13_0 - v.r13_1 - v.r13_2 - v.r13_3 - v.r13_4 - v.r13_5 - v.r13_6 + v.r29;
    dStates.Ops = v.r6_1 + v.r6_2 + v.r6_3 + v.r6_4 + v.r6_5 + v.r6_6 + v.r7_0 + v.r7_1 + v.r7_2 + v.r7_3 + v.r7_4 + v.r7_5 + v.r7_6 - v.r8 + v.r11 - v.r12;
    dStates.Ops_G = v.r9 - v.r10;
    dStates.Ops_G_GTP = v.r10 - v.r11;
    dStates.Ops_Gt = v.r8 - v.r9;
    dStates.PDE = -v.r18 + v.r25 + v.r26;
    dStates.PDE_Ga_GTP = v.r18 - v.r19;
    dStates.PDE_a_Ga_GTP = v.r19 - v.r20 + v.r23 - v.r24 - v.r26 + v.r27;
    dStates.R = -v.r1 + v.r12 - v.GtRpre;
    dStates.R0 = v.r1 - v.r2_0 - v.r7_0 - v.r13_0 + v.r16_0;
    dStates.R0_G = v.r14_0 - v.r15_0;
    dStates.R0_G_GTP = v.r15_0 - v.r16_0;
    dStates.R0_Gt = v.rstprec + v.r13_0 - v.r14_0;
    dStates.R0_RKpre = v.r2_0 - v.r3_0;
    dStates.R1 = -v.r2_1 + v.r4_1 - v.r5_1 - v.r7_1 - v.r13_1 + v.r16_1;
    dStates.R1_Arr = v.r5_1 - v.r6_1;
    dStates.R1_G = v.r14_1 - v.r15_1;
    dStates.R1_G_GTP = v.r15_1 - v.r16_1;
    dStates.R1_Gt = v.r13_1 - v.r14_1;
    dStates.R1_RKpost = v.r3_0 - v.r4_1;
    dStates.R1_RKpre = v.r2_1 - v.r3_1;
    dStates.R2 = -v.r2_2 + v.r4_2 - v.r5_2 - v.r7_2 - v.r13_2 + v.r16_2;
    dStates.R2_Arr = v.r5_2 - v.r6_2;
    dStates.R2_G = v.r14_2 - v.r15_2;
    dStates.R2_G_GTP = v.r15_2 - v.r16_2;
    dStates.R2_Gt = v.r13_2 - v.r14_2;
    dStates.R2_RKpost = v.r3_1 - v.r4_2;
    dStates.R2_RKpre = v.r2_2 - v.r3_2;
    dStates.R3 = -v.r2_3 + v.r4_3 - v.r5_3 - v.r7_3 - v.r13_3 + v.r16_3;
    dStates.R3_Arr = v.r5_3 - v.r6_3;
    dStates.R3_G = v.r14_3 - v.r15_3;
    dStates.R3_G_GTP = v.r15_3 - v.r16_3;
    dStates.R3_Gt = v.r13_3 - v.r14_3;
    dStates.R3_RKpost = v.r3_2 - v.r4_3;
    dStates.R3_RKpre = v.r2_3 - v.r3_3;
    dStates.R4 = -v.r2_4 + v.r4_4 - v.r5_4 - v.r7_4 - v.r13_4 + v.r16_4;
    dStates.R4_Arr = v.r5_4 - v.r6_4;
    dStates.R4_G = v.r14_4 - v.r15_4;
    dStates.R4_G_GTP = v.r15_4 - v.r16_4;
    dStates.R4_Gt = v.r13_4 - v.r14_4;
    dStates.R4_RKpost = v.r3_3 - v.r4_4;
    dStates.R4_RKpre = v.r2_4 - v.r3_4;
    dStates.R5 = -v.r2_5 + v.r4_5 - v.r5_5 - v.r7_5 - v.r13_5 + v.r16_5;
    dStates.R5_Arr = v.r5_5 - v.r6_5;
    dStates.R5_G = v.r14_5 - v.r15_5;
    dStates.R5_G_GTP = v.r15_5 - v.r16_5;
    dStates.R5_Gt = v.r13_5 - v.r14_5;
    dStates.R5_RKpost = v.r3_4 - v.r4_5;
    dStates.R5_RKpre = v.r2_5 - v.r3_5;
    dStates.R6 = -v.r2_6 + v.r4_6 - v.r5_6 - v.r7_6 - v.r13_6 + v.r16_6;
    dStates.R6_Arr = v.r5_6 - v.r6_6;
    dStates.R6_G = v.r14_6 - v.r15_6;
    dStates.R6_G_GTP = v.r15_6 - v.r16_6;
    dStates.R6_Gt = v.r13_6 - v.r14_6;
    dStates.R6_RKpost = v.r3_5 - v.r4_6;
    dStates.R6_RKpre = v.r2_6;
    dStates.RGS = -v.r22 + v.r23 - v.r24 + v.r25;
    dStates.RGS_Ga_GTP_a_PDE_a_Ga_GTP = v.r22 - v.r23;
    dStates.RGS_PDE_a_Ga_GTP = v.r24 - v.r25;
    dStates.RK = -v.r2_0 - v.r2_1 - v.r2_2 - v.r2_3 - v.r2_4 - v.r2_5 - v.r2_6 + v.r4_1 + v.r4_2 + v.r4_3 + v.r4_4 + v.r4_5 + v.r4_6 - v.r31;
    dStates.R_Gt = -v.rstprec + v.GtRpre;
    dStates.RecR_Ca = v.r30 - v.r31;
    dStates.RecR_Ca_RK = v.r31;
    dStates.RecT = -v.r30;
    
    %% 5. 计算ATP消耗
    v.ATPtrans = v.r3_0 + v.r3_1 + v.r3_2 + v.r3_3 + v.r3_4 + v.r3_5 + ...
                 v.r10 + v.r15_0 + v.r15_1 + v.r15_2 + v.r15_3 + v.r15_4 + ...
                 v.r15_5 + v.r15_6;
end