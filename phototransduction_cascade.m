% phototransduction_cascade.m
function [dStates, E_star, v] = phototransduction_cascade(states, stimulus, params)
% Calculates the derivatives (d/dt) of all state variables in the phototransduction cascade.
%
% Outputs:
%   dStates - A struct containing the derivatives (d/dt) of all state variables.
%   E_star  - The concentration of activated PDE.
%   v       - A struct containing the rates of all reactions.

    %% 1. Calculate Rate Constants
    % Phosphorylation-dependent rate constants
    kRK1_1 = params.kRK1_0 * exp(-params.omega); kRK1_2 = params.kRK1_0 * exp(-params.omega * 2);
    kRK1_3 = params.kRK1_0 * exp(-params.omega * 3); kRK1_4 = params.kRK1_0 * exp(-params.omega * 4);
    kRK1_5 = params.kRK1_0 * exp(-params.omega * 5); kRK1_6 = 0;
    
    % Arrestin binding rate constants
    kA1_1 = params.kArr; kA1_2 = params.kArr + 1 * params.m_Arr; kA1_3 = params.kArr + 2 * params.m_Arr;
    kA1_4 = params.kArr + 3 * params.m_Arr; kA1_5 = params.kArr + 3 * params.m_Arr; kA1_6 = params.kArr + 3 * params.m_Arr;
    
    % G-protein binding rate constants
    kGpre1 = params.kG1_0 * 1.6; kGpre2 = params.kG2 * 315;
    kG1_1 = params.kG1_0 * exp(-params.omega_G); kG1_2 = params.kG1_0 * exp(-params.omega_G * 2);
    kG1_3 = params.kG1_0 * exp(-params.omega_G * 3); kG1_4 = params.kG1_0 * exp(-params.omega_G * 4);
    kG1_5 = params.kG1_0 * exp(-params.omega_G * 5); kG1_6 = params.kG1_0 * exp(-params.omega_G * 6);
    
    %% 2. Calculate Reaction Rates (v)
    v.r1 = stimulus * states.R / params.Rtot; % Rhodopsin activation by light
    v.rstprec = stimulus * states.R_Gt / params.Rtot;

    % R-RK Reactions (Rhodopsin Kinase)
    v.r2_0 = params.kRK1_0 * states.RK * states.R0 - params.kRK2 * states.R0_RKpre;
    v.r2_1 = kRK1_1 * states.RK * states.R1 - params.kRK2 * states.R1_RKpre;
    % ... (and so on for R2-R6)
    v.r2_2 = kRK1_2 * states.RK * states.R2 - params.kRK2 * states.R2_RKpre;
    v.r2_3 = kRK1_3 * states.RK * states.R3 - params.kRK2 * states.R3_RKpre;
    v.r2_4 = kRK1_4 * states.RK * states.R4 - params.kRK2 * states.R4_RKpre;
    v.r2_5 = kRK1_5 * states.RK * states.R5 - params.kRK2 * states.R5_RKpre;
    v.r2_6 = kRK1_6 * states.RK * states.R6 - params.kRK2 * states.R6_RKpre;

    % Phosphorylation steps
    v.r3_0 = params.kRK3_ATP * states.R0_RKpre;
    v.r3_1 = params.kRK3_ATP * states.R1_RKpre;
    % ... (and so on for R2-R5)
    v.r3_2 = params.kRK3_ATP * states.R2_RKpre; v.r3_3 = params.kRK3_ATP * states.R3_RKpre;
    v.r3_4 = params.kRK3_ATP * states.R4_RKpre; v.r3_5 = params.kRK3_ATP * states.R5_RKpre;
    
    % RK dissociation
    v.r4_1 = params.kRK4 * states.R1_RKpost;
    v.r4_2 = params.kRK4 * states.R2_RKpost;
    % ... (and so on for R3-R6)
    v.r4_3 = params.kRK4 * states.R3_RKpost; v.r4_4 = params.kRK4 * states.R4_RKpost;
    v.r4_5 = params.kRK4 * states.R5_RKpost; v.r4_6 = params.kRK4 * states.R6_RKpost;

    % Arrestin Reactions
    v.r5_1 = kA1_1 * states.Arr * states.R1 - params.kA2 * states.R1_Arr;
    v.r5_2 = kA1_2 * states.Arr * states.R2 - params.kA2 * states.R2_Arr;
    % ... (and so on for R3-R6)
    v.r5_3 = kA1_3 * states.Arr * states.R3 - params.kA2 * states.R3_Arr; v.r5_4 = kA1_4 * states.Arr * states.R4 - params.kA2 * states.R4_Arr;
    v.r5_5 = kA1_5 * states.Arr * states.R5 - params.kA2 * states.R5_Arr; v.r5_6 = kA1_6 * states.Arr * states.R6 - params.kA2 * states.R6_Arr;

    v.r6_1 = params.kA3 * states.R1_Arr; v.r6_2 = params.kA3 * states.R2_Arr;
    % ... (and so on for R3-R6)
    v.r6_3 = params.kA3 * states.R3_Arr; v.r6_4 = params.kA3 * states.R4_Arr;
    v.r6_5 = params.kA3 * states.R5_Arr; v.r6_6 = params.kA3 * states.R6_Arr;
    
    % Thermal Activation
    v.r7_0 = params.ktherm * states.R0; v.r7_1 = params.ktherm * states.R1;
    % ... (and so on for R2-R6)
    v.r7_2 = params.ktherm * states.R2; v.r7_3 = params.ktherm * states.R3; v.r7_4 = params.ktherm * states.R4;
    v.r7_5 = params.ktherm * states.R5; v.r7_6 = params.ktherm * states.R6;

    % G-protein (Transducin) Reactions
    v.r8 = params.kOps * states.Ops * states.Gt - params.kG2 * states.Ops_Gt;
    v.r9 = params.kG3 * states.Ops_Gt - params.kG4_GDP * states.Ops_G;
    v.r10 = params.kG5_GTP * states.Ops_G; v.r11 = params.kG6 * states.Ops_G_GTP;
    v.r12 = params.kRrecyc * states.Ops; % Opsin recycling
    
    v.GtRpre = kGpre1 * states.Gt * states.R - kGpre2 * states.R_Gt;
    
    v.r13_0 = params.kG1_0 * states.Gt * states.R0 - params.kG2 * states.R0_Gt;
    % ... (and so on for R1-R6)
    v.r13_1 = kG1_1 * states.Gt * states.R1 - params.kG2 * states.R1_Gt; v.r13_2 = kG1_2 * states.Gt * states.R2 - params.kG2 * states.R2_Gt;
    v.r13_3 = kG1_3 * states.Gt * states.R3 - params.kG2 * states.R3_Gt; v.r13_4 = kG1_4 * states.Gt * states.R4 - params.kG2 * states.R4_Gt;
    v.r13_5 = kG1_5 * states.Gt * states.R5 - params.kG2 * states.R5_Gt; v.r13_6 = kG1_6 * states.Gt * states.R6 - params.kG2 * states.R6_Gt;

    v.r14_0 = params.kG3 * states.R0_Gt - params.kG4_GDP * states.R0_G;
    % ... (and so on for R1-R6)
    v.r14_1 = params.kG3 * states.R1_Gt - params.kG4_GDP * states.R1_G; v.r14_2 = params.kG3 * states.R2_Gt - params.kG4_GDP * states.R2_G;
    v.r14_3 = params.kG3 * states.R3_Gt - params.kG4_GDP * states.R3_G; v.r14_4 = params.kG3 * states.R4_Gt - params.kG4_GDP * states.R4_G;
    v.r14_5 = params.kG3 * states.R5_Gt - params.kG4_GDP * states.R5_G; v.r14_6 = params.kG3 * states.R6_Gt - params.kG4_GDP * states.R6_G;

    v.r15_0 = params.kG5_GTP * states.R0_G;
    % ... (and so on for R1-R6)
    v.r15_1 = params.kG5_GTP * states.R1_G; v.r15_2 = params.kG5_GTP * states.R2_G; v.r15_3 = params.kG5_GTP * states.R3_G;
    v.r15_4 = params.kG5_GTP * states.R4_G; v.r15_5 = params.kG5_GTP * states.R5_G; v.r15_6 = params.kG5_GTP * states.R6_G;

    v.r16_0 = params.kG6 * states.R0_G_GTP;
    % ... (and so on for R1-R6)
    v.r16_1 = params.kG6 * states.R1_G_GTP; v.r16_2 = params.kG6 * states.R2_G_GTP; v.r16_3 = params.kG6 * states.R3_G_GTP;
    v.r16_4 = params.kG6 * states.R4_G_GTP; v.r16_5 = params.kG6 * states.R5_G_GTP; v.r16_6 = params.kG6 * states.R6_G_GTP;

    v.r17 = params.kG7 * states.G_GTP; % G-alpha-GTP production
    
    % PDE Reactions
    v.r18 = params.kP1 * states.PDE * states.Ga_GTP - params.kP1_rev * states.PDE_Ga_GTP;
    v.r19 = params.kP2 * states.PDE_Ga_GTP;
    v.r20 = params.kP3 * states.PDE_a_Ga_GTP * states.Ga_GTP;
    v.r21 = params.kP4 * states.Ga_GTP_PDE_a_Ga_GTP;
    
    % RGS Reactions
    v.r22 = params.kRGS1 * states.RGS * states.Ga_GTP_a_PDE_a_Ga_GTP;
    v.r23 = params.kRGS2 * states.RGS_Ga_GTP_a_PDE_a_Ga_GTP;
    v.r24 = params.kRGS1 * states.RGS * states.PDE_a_Ga_GTP;
    v.r25 = params.kRGS2 * states.RGS_PDE_a_Ga_GTP;
    
    % Shutoff and Recycling
    v.r26 = params.kPDEshutoff * states.PDE_a_Ga_GTP;
    v.r27 = params.kPDEshutoff * states.Ga_GTP_a_PDE_a_Ga_GTP;
    v.r28 = params.kGshutoff * states.Ga_GTP;
    v.r29 = params.kGrecyc * states.Gbg * states.Ga_GDP; % G-protein recycling
    
    % Recoverin Reactions
    v.r30 = params.kRec1 * states.RecT * states.Caos - params.kRec2 * states.RecR_Ca;
    v.r31 = params.kRec3 * states.RecR_Ca * states.RK - params.kRec4 * states.RecR_Ca_RK;
    
    % Arrestin Dimerization
    v.r_diarr = params.kA4 * states.Arr * states.Arr - params.kA5 * states.Arr_di;
    v.r_tetraarr = params.kA4 * states.Arr_di * states.Arr_di - params.kA5 * states.Arr_tetra;
    
    %% 3. Calculate E_star (Activated PDE)
    E_star = states.PDE_a_Ga_GTP + 2 * states.Ga_GTP_a_PDE_a_Ga_GTP + states.Ga_GTP_PDE_a_Ga_GTP;
    
    %% 4. Calculate State Derivatives (d/dt)
    dStates = struct();
    % The following equations define the rate of change for each species in the cascade.
    dStates.Arr = (-v.r5_1 - v.r5_2 - v.r5_3 - v.r5_4 - v.r5_5 - v.r5_6 + v.r6_1 + v.r6_2 + v.r6_3 + v.r6_4 + v.r6_5 + v.r6_6 - 2*v.r_diarr);
    dStates.Arr_di = v.r_diarr - 2*v.r_tetraarr;
    dStates.Arr_tetra = v.r_tetraarr;
    dStates.G_GTP = v.r11 + sum([v.r16_0, v.r16_1, v.r16_2, v.r16_3, v.r16_4, v.r16_5, v.r16_6]) - v.r17;
    dStates.Ga_GDP = v.r23 + v.r25 + v.r26 + v.r27 + v.r28 - v.r29;
    dStates.Ga_GTP = v.r17 - v.r18 - v.r20 - v.r28;
    dStates.Ga_GTP_PDE_a_Ga_GTP = v.r20 - v.r21;
    dStates.Ga_GTP_a_PDE_a_Ga_GTP = v.r21 - v.r22 - v.r27;
    dStates.Gbg = v.r17 - v.r29;
    dStates.Gt = -v.r8 - v.GtRpre - sum([v.r13_0, v.r13_1, v.r13_2, v.r13_3, v.r13_4, v.r13_5, v.r13_6]) + v.r29;
    dStates.Ops = sum([v.r6_1, v.r6_2, v.r6_3, v.r6_4, v.r6_5, v.r6_6]) + sum([v.r7_0, v.r7_1, v.r7_2, v.r7_3, v.r7_4, v.r7_5, v.r7_6]) - v.r8 + v.r11 - v.r12;
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
    % ... (and so on for all other R states and complexes)
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
    dStates.RK = -sum([v.r2_0, v.r2_1, v.r2_2, v.r2_3, v.r2_4, v.r2_5, v.r2_6]) + sum([v.r4_1, v.r4_2, v.r4_3, v.r4_4, v.r4_5, v.r4_6]) - v.r31;
    dStates.R_Gt = -v.rstprec + v.GtRpre;
    dStates.RecR_Ca = v.r30 - v.r31;
    dStates.RecR_Ca_RK = v.r31;
    dStates.RecT = -v.r30;
    
    %% 5. Calculate ATP Consumption Rate for this cascade
    v.ATPtrans = sum([v.r3_0, v.r3_1, v.r3_2, v.r3_3, v.r3_4, v.r3_5, v.r10]) + ...
                 sum([v.r15_0, v.r15_1, v.r15_2, v.r15_3, v.r15_4, v.r15_5, v.r15_6]);
end