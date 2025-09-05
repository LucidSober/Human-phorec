function states_r = Y2States(Y)

    % 电生理状态 (Gating)
    states_r.Vm = Y(1);
    states_r.p_CaL = Y(2);
    states_r.p_h = Y(3);
    states_r.p_Kv = Y(4);
    states_r.m_KCa = Y(5);

    % 离子浓度状态 (μM for Ca, mM for others)
    states_r.Ki = Y(6);
    states_r.Nai = Y(7);
    states_r.Cli = Y(8);
    states_r.Caos = Y(9);
    states_r.Cais = Y(10);
    states_r.Caif = Y(11);

    % 钙离子缓冲状态 (μM)
    states_r.Ca_buff_os = Y(12);
    states_r.Ca_ls = Y(13);
    states_r.Ca_hs = Y(14);
    states_r.Ca_lf = Y(15);
    states_r.Ca_hf = Y(16);
    states_r.Cab = Y(17);

    % 光转导状态
    states_r.cGMP = Y(18);
    states_r.PDE_active = Y(19); % 遗留占位符
    states_r.Arr = Y(20); states_r.Arr_di = Y(21); states_r.Arr_tetra = Y(22);
    states_r.G_GTP = Y(23); states_r.Ga_GDP = Y(24); states_r.Ga_GTP = Y(25);
    states_r.Ga_GTP_PDE_a_Ga_GTP = Y(26); states_r.Ga_GTP_a_PDE_a_Ga_GTP = Y(27);
    states_r.Gbg = Y(28); states_r.Gt = Y(29);
    states_r.Ops = Y(30); states_r.Ops_G = Y(31); states_r.Ops_G_GTP = Y(32); states_r.Ops_Gt = Y(33);
    states_r.PDE = Y(34); states_r.PDE_Ga_GTP = Y(35); states_r.PDE_a_Ga_GTP = Y(36);
    states_r.R = Y(37); states_r.R0 = Y(38); states_r.R0_G = Y(39); states_r.R0_G_GTP = Y(40);
    states_r.R0_Gt = Y(41); states_r.R0_RKpre = Y(42);
    states_r.R1 = Y(43); states_r.R1_Arr = Y(44); states_r.R1_G = Y(45); states_r.R1_G_GTP = Y(46);
    states_r.R1_Gt = Y(47); states_r.R1_RKpost = Y(48); states_r.R1_RKpre = Y(49);
    states_r.R2 = Y(50); states_r.R2_Arr = Y(51); states_r.R2_G = Y(52); states_r.R2_G_GTP = Y(53);
    states_r.R2_Gt = Y(54); states_r.R2_RKpost = Y(55); states_r.R2_RKpre = Y(56);
    states_r.R3 = Y(57); states_r.R3_Arr = Y(58); states_r.R3_G = Y(59); states_r.R3_G_GTP = Y(60);
    states_r.R3_Gt = Y(61); states_r.R3_RKpost = Y(62); states_r.R3_RKpre = Y(63);
    states_r.R4 = Y(64); states_r.R4_Arr = Y(65); states_r.R4_G = Y(66); states_r.R4_G_GTP = Y(67);
    states_r.R4_Gt = Y(68); states_r.R4_RKpost = Y(69); states_r.R4_RKpre = Y(70);
    states_r.R5 = Y(71); states_r.R5_Arr = Y(72); states_r.R5_G = Y(73); states_r.R5_G_GTP = Y(74);
    states_r.R5_Gt = Y(75); states_r.R5_RKpost = Y(76); states_r.R5_RKpre = Y(77);
    states_r.R6 = Y(78); states_r.R6_Arr = Y(79); states_r.R6_G = Y(80); states_r.R6_G_GTP = Y(81);
    states_r.R6_Gt = Y(82); states_r.R6_RKpost = Y(83); states_r.R6_RKpre = Y(84);
    states_r.RGS = Y(85); states_r.RGS_Ga_GTP_a_PDE_a_Ga_GTP = Y(86); states_r.RGS_PDE_a_Ga_GTP = Y(87);
    states_r.RK = Y(88); states_r.R_Gt = Y(89);
    states_r.RecT = Y(90); states_r.RecR_Ca = Y(91); states_r.RecR_Ca_RK = Y(92);

    states_r.Rtot = Y(93);
    states_r.E_star = Y(94);

end