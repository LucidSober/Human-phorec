%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function I = corrected_RC_stim(t)
% 修正的双段RC刺激函数
% 第一段: 20μA + 20μA*(1-exp(...)) → 40μA (1.25s到1.75s)
% 第二段: 20μA*exp(...) → 0μA (1.75s之后)
% 
% 输入:
%   t - 时间向量 (s)
% 输出:
%   I - 电流 (A)

% 参数定义
I_offset = 2e-5;      % 偏移量 20 μA
I_amp = 2e-5;         % 幅度 20 μA
t_start = 1.25;       % 起始时间 (s)
t_mid = 1.75;         % 衔接时间 (s)
tau1 = 0.167;         % 第一段时间常数 (s)
tau2 = 0.361;         % 第二段时间常数 (s)

% 初始化输出
I = zeros(size(t));

% 第一段：偏移 + 反向RC充电 (1.25s ≤ t ≤ 1.75s)
idx1 = (t >= t_start) & (t <= t_mid);
I(idx1) = I_offset + I_amp * (1 - exp(-(t_mid - t(idx1)) / tau1));

% 第二段：RC放电 (t > 1.75s)
idx2 = t > t_mid;
I(idx2) = I_offset * exp(-(t(idx2) - t_mid) / tau2);

end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
