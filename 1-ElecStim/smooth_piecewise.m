%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function y = smooth_piecewise(te_real,start_real,ms,me,ys,ye)
    global params
    
    % start=params.t_amp_top_real;
    L  = params.t_Herm_dt;
    % amp=params.stim_amp;
    
    m0 = ms;  m1 = me;    % 端点斜率 (对 x)
    y0 = ys;  y1 = ye; % 端点值
    
    s = (te_real - start_real) / L;

    h00 = 2*s.^3 - 3*s.^2 + 1;
    h10 = s.^3   - 2*s.^2 + s;
    h01 = -2*s.^3 + 3*s.^2;
    h11 = s.^3   - s.^2;

    y = h00.*y0 + h10.*(m0*L) + h01.*y1 + h11.*(m1*L);


end

%[appendix]{"version":"1.0"}
%---
