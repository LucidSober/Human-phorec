%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function I = corrected_RC_stim(t)
% 修正的双段RC刺激函数
% 第一段: 20μA + 20μA*(1-exp(...)) → 40μA (1.25s到1.75s)
% 第二段: 20μA*exp(...) → 0μA (1.75s之后)
% 
% 输入:
%   t -  (s)
% 输出:
%   I - 电流 (A)

global params


% 参数定义
I_offset = params.stim_amp/2;      % 偏移量 
I_amp = params.stim_amp/2;         % 幅度 
t_start = params.t_amp_top_real;       % 起始时间 (s)
t_mid = params.t_amp_top_real+0.5;         % 衔接时间 (s)@@@
tau1 = 0.167;         % 第一段时间常数 (s)
tau2 = 0.361;         % 第二段时间常数 (s)


% 第一段：偏移 + 反向RC充电 (1.25s ≤ t ≤ 1.75s)
% idx1 = (t >= t_start) & (t <= t_mid);

% if( (t >= t_start) && (t <= t_mid) )
%     I = I_offset + I_amp * (1 - exp(-(t_mid - t) / tau1));
% else
%     I = I_offset * exp( (t_mid-t)/tau2 );
% end

I = I_offset + I_amp * (1 - exp(-(t_mid - t) / tau1));

% I = I_offset/2 * exp( (t_mid-t)/tau2 );




end

%[appendix]{"version":"1.0"}
%---
