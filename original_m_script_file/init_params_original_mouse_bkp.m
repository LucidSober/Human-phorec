function params = init_params()
    % 返回包含所有参数的结构体
    
    %% 光转导级联参数
    % RK (rhodopsin kinase) 参数
    params.kRK1_0 = 0.1724;
    params.omega = 2.5;
    params.kRK2 = 250;
    params.kRK3_ATP = 4000;
    params.kRK4 = 250;
    
    % Arrestin 参数
    params.kArr = 0.0000099147;
    params.kA2 = 0.026;
    params.m_Arr = 0.0000095475;
    params.kA3 = 1.1651;
    params.kA4 = 0.00000029965;
    params.kA5 = 0.424;
    
    % Opsin 参数
    params.kOps = 0.00000000000061172;
    params.kRrecyc = 0.0007;
    
    % G蛋白参数
    params.omega_G = 0.6;
    params.kG1_0 = 0.001;
    params.kG2 = 2200;
    params.kG3 = 8500;
    params.kG4_GDP = 400;
    params.kG5_GTP = 3500;
    params.kG6 = 8500;
    params.kG7 = 200;
    params.kGrecyc = 2;
    params.kGshutoff = 0.05;
    
    % PDE参数
    params.kP1 = 0.05497;
    params.kP1_rev = 0;
    params.kP2 = 940.7;
    params.kP3 = 0.0000000014983;
    params.kP4 = 21.088;
    params.kPDEshutoff = 0.1;
    
    % RGS参数
    params.kRGS1 = 0.000048182;
    params.kRGS2 = 98;
    
    % Recoverin参数
    params.kRec1 = 0.011;
    params.kRec2 = 0.05;
    params.kRec3 = 0.00041081;
    params.kRec4 = 0.610084;
    
    % 热激活
    params.ktherm = 0.0238;
    
    %% cGMP和钙参数
    params.Vcyto = 0.03916; % pL
    params.Kc1 = 0.171; % μM
    params.Kc2 = 0.059; % μM
    params.m1 = 3;
    params.m2 = 1.5;
    params.alfamax = 60; % μM/s
    params.betadark = 3.1872748181298; % s^-1
    params.betasub = 0.021826; % s^-1/molecule
    
    % 钙相关参数
    params.fCa = 0.12;
    params.Jdark = 14.87; % pA
    params.cGMPdark = 6.4944; % μM
    params.ncg = 3.8;
    params.gammaCa = 981.3558; % s^-1
    params.Caos_dark = 0.25; % μM
    params.Caos_0 = 0.023; % μM
    
    % 钙缓冲参数
    params.k1 = 9.37059; % μM^-1 s^-1
    params.k2 = 46.412; % s^-1
    params.eT = 400; % μM
    params.k3 = 1e-20; % s^-1 (钙扩散)
    
    %% 细胞参数
    params.Cm = 3.6; % pF
    params.F = 96485.34; % C/mol
    params.R_gas = 8.3145; % J/(mol·K)
    params.T = 310; % K
    params.NA = 6.022e23; % molecules/mol
    
    % 体积参数
    params.V_os = 0.00000000000003916; % L (外段)
    params.V_is = 0.000000000000023496; % L (内段)
    params.V_if = 0.000000000000015664; % L (内段纤维)
    params.V_cell = params.V_os + params.V_is + params.V_if;
    params.Vos = params.V_os;
    params.Vis = params.V_is;
    params.Vif = params.V_if;
    
    %% 离子浓度
    % 外部浓度
    params.Cao = 1.6; % mM (1600 μM = 1.6 mM)
    params.Nao = 145; % mM
    params.Ko = 5; % mM
    params.Clo = 110; % mM
    
    % 初始内部浓度（用于参考）
    params.Nai_0 = 5; % mM
    params.Ki_0 = 140; % mM
    params.Cli_0 = 30; % mM
    params.Cais_0 = 0.1; % μM
    
    %% 钙系统参数
    params.DCa = 6e-8; % dm^2/s
    params.delta = 3e-5; % dm
    params.S_1 = 3.142e-8; % dm^2
    
    % 钙缓冲池参数
    params.Lb1 = 0.4; % s^-1 μM^-1
    params.Lb2 = 0.21; % s^-1
    params.Hb1 = 100; % s^-1 μM^-1
    params.Hb2 = 90; % s^-1
    params.BL = 500; % μM
    params.BH = 300; % μM
    
    %% 通道电导
    params.gKv = 0.2; % nS
    params.gCl = 2; % nS
    params.gKCa = 0.113; % nS

    % Ih 通道通透性 (来自补充材料 S3)
    params.P_K_h = 0.15;
    params.P_Na_h = 0.05;
    
    % ICaL 通道通透性 (来自补充材料 S3)
    % 首先定义一个缩放因子
    k_CaL_scale_factor = 0.65;
    params.P_Ca_CaL = 1 * k_CaL_scale_factor;
    params.P_K_CaL = 0.000365 * k_CaL_scale_factor;
    params.P_Na_CaL = 0.0000185 * k_CaL_scale_factor;
    
    % CNG 通道通透性 (来自补充材料 S3)
    params.P_K_CNG = 1 / 3.6;
    params.P_Na_CNG = 1 / 3.6;
    params.P_Ca_CNG = 6 / 3.6;
    
    %% 泵和交换器参数
    params.kPMCA = 10; % pA/μM
    params.Km_Cai = 0.0005; % mM PMCA的钙亲和度
    params.kNCX = 250; % pA
    params.k_NKCC1 = 0.005; % s^-1
    params.k_KCC2 = 0.046; % s^-1
    params.kNKCC1 = params.k_NKCC1;
    params.kKCC2 = params.k_KCC2;
    
    %% 泄漏电导
    params.Gleak_K = 0.00011488937611594242; % nS
    params.Gleak_Na = 0.000071150566012260358; % nS
    params.Gleak_Cais = 0.0002481452439628276; % nS
    params.Gleak_Cl = 0.00035986965484864007; % nS
    params.Gleak_Caos = 0.000021232320024139139; % nS
    
    %% 数值计算参数
    params.numConcFactor = 1 / (6.022e5 * params.Vcyto); % 浓度转换因子
    
    %% 总量常数（用于光转导）
    params.Rtot = 100000000.0;
    params.PDEtot = 2000000.0;
    params.Gtot = 10000000.0;
    params.RGStot = 100000.0;
    params.ArrTot = 7074600.0;
end