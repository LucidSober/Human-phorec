%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function [y, k] = smooth_piecewise_back(te_abs, x1, x2, y1, beta_up, beta_dn)
% 起点 (x1,y1), 初始斜率=0
% 先缓慢下降到恒定斜率 k(<0)，中段以 k 直线下降，
% 末段在到达 (x2,0) 前斜率缓慢回到 0，终点斜率=0。
%
% te_abs 可为标量或向量，范围 [x1,x2]

    if nargin < 5, beta_up = 0.3; end
    if nargin < 6, beta_dn = 0.2; end

    L = x2 - x1;
    if ~(L > 0), error('需要 x2 > x1'); end

    % 归一化 s ∈ [0,1]
    s  = (te_abs - x1) / L;
    sc = max(0, min(1, s));

    % 夹紧并避免退化
    b1 = min(max(beta_up, 1e-9), 1-1e-9);
    b2 = min(max(beta_dn, 1e-9), 1-1e-9);
    if b1 + b2 >= 1
        scale = 0.999 / (b1 + b2);
        b1 = b1*scale; b2 = b2*scale;
    end

    % ---------- 构造无量纲积分 F(s) ----------
    % 斜率 m(s) = k * g(s)
    % 段A 0..b1:       g = smoothstep(s/b1)
    % 段B b1..1-b2:    g = 1
    % 段C 1-b2..1:     g = 1 - smoothstep((s-(1-b2))/b2)
    F = zeros(size(sc));

    % 段A: ∫ smoothstep = s^3/b1^2 - s^4/(2 b1^3)
    maskA = sc <= b1;
    sa = sc(maskA);
    Fa = sa.^3/(b1^2) - sa.^4/(2*b1^3);

    % 段B: 承接 F(b1)=b1/2，之后线性增长
    maskB = (sc > b1) & (sc <= 1-b2);
    sb = sc(maskB);
    Fb = (sb - b1/2);

    % 段C: 承接 F(1-b2)=1-b2 - b1/2，再加 ∫(1 - smoothstep)
    maskC = sc > 1-b2;
    scC = sc(maskC);
    w = (scC - (1-b2)) / b2;                         % w∈(0,1]
    Fc = (1 - b2 - b1/2) + b2 .* ( w - w.^3 + 0.5*w.^4 );

    F(maskA) = Fa;
    F(maskB) = Fb;
    F(maskC) = Fc;

    % 终点条件 y(x2)=0 ⇒ y1 + L*k*F(1) = 0
    F1 = 1 - (b1 + b2)/2;

    % —— 强制 k 为负值 ——（无论 y1 正负，均取负号）
    k = - abs( y1 / (L * F1) );

    % 位置
    y = y1 + L * k .* F;
end




% function [y, k] = smooth_piecewise_back(te_real, x1, x2, y1, beta)
% % 起点 (x1,y1), 初始斜率 k0=0
% % 终点 (x2,0), 斜率为未知 k（由边界条件自动求解）
% % 先把斜率从 0 用 smoothstep 平滑提升到 k，再以常斜率 k 线性到终点
% %
% % 输入:
% %   te_abs : 标量或向量，x ∈ [x1,x2]
% %   x1, x2 : 起止横坐标 (x2 > x1)
% %   y1     : 起点纵坐标
% %   beta   : 过渡比例 (0~1)，建议 0.2~0.5；可省略，默认 0.3
% %
% % 输出:
% %   y : y(te_abs)
% %   k : 过渡后的稳定斜率（负值表示下降）
% 
%     if nargin < 5, beta = 0.3; end
%     L = x2 - x1;
%     if ~(L > 0), error('需要 x2 > x1'); end
% 
%     % 归一化到 s∈[0,1]
%     s  = (te_real - x1) / L;
%     sc = max(0, min(1, s));
% 
%     % 防止 beta=0 或 1
%     b = min(max(beta, 1e-9), 1 - 1e-9);
% 
%     % ---- 斜率轨迹 m(s) 的积分 F(s)（使 y = y1 + L*k*F(s)）----
%     % 0<=s<=b: m(s)=k*smoothstep(s/b),  smoothstep(u)=3u^2-2u^3
%     % 积分得：F(s)= s^3/b^2 - s^4/(2 b^3)
%     % s>b:    m(s)=k,          F(s)= s - b/2
%     F = zeros(size(sc));
%     mask = sc <= b;
%     sb = sc(mask);
%     F(mask)  = sb.^3/(b^2) - sb.^4/(2*b^3);
%     F(~mask) = sc(~mask) - b/2;
% 
%     % ---- 用终点 y(x2)=0 解 k ----
%     % F(1) = 1 - b/2
%     k = - y1 / (L * (1 - b/2));
% 
%     % ---- 计算 y(t) ----
%     y = y1 + L * k .* F;
% end



%[appendix]{"version":"1.0"}
%---
