%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function I = I_func(te)

global params

f = params.f; 
amp=params.stim_amp;
pt=params.stim_keypoint_t;
lp=params.stim_linparam;

pha = mod(te, 1/f);

if( pha>=pt(1) && pha<=pt(2) )
    I = lp(1,1) * pha + lp(1,2);
elseif( pha>pt(2) && pha<pt(3) )
    I = amp;
elseif( pha>=pt(3) && pha<=pt(4) )
    I = lp(2,1) * pha + lp(2,2);
elseif( pha>pt(4) && pha<pt(5) )
    I = 0;
elseif( pha>=pt(5) && pha<=pt(6) )
    I = lp(3,1) * pha + lp(3,2);
elseif( pha>pt(6) && pha<pt(7) )
    I = -amp;
elseif( pha>=pt(7) && pha<=pt(8) )
    I = lp(4,1) * pha + lp(4,2);
else
    I = 0;
end



% I=interp1(params.min_step_t, params.Stim_wave, te, 'previous', 0);

end

%[appendix]{"version":"1.0"}
%---
