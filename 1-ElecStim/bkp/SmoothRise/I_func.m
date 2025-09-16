%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function I = I_func(te)

global params 

f = params.f; 
amp=params.stim_amp;
pt=params.stim_keypoint_t;
lp=params.stim_linparam;

pha = mod(te, 1/f);

if(pha<0)  
   I = 0;    
elseif(pha < (pt(2)-pt(1)))
    I = smooth_piecewise(pha+(params.stim_start+pt(1)), params.stim_start+pt(1), 0,params.k0_Herm,0,params.y0_Herm); 
    % I = 0;
elseif( pha >= (pt(2)-pt(1)) )
    pha=pha+pt(1);

    if( pha>=pt(2) && pha<pt(3) )
        I = lp(1,1) * pha + lp(1,2);
        % I = 0;
    elseif( pha>=pt(3) && pha<=pt(4) )
        I = smooth_piecewise(pha+params.stim_start, params.t_amp_top_real, params.k0_Herm,0,amp,params.y1_Herm);
        % I =0;
    elseif( pha>pt(4) && pha<=pt(5) )
        I = lp(2,1) * pha + lp(2,2);
        % I =0;
    else
        I = 0;
    end

end


% if( pha>=pt(1) && pha<pt(2) )
%     I = lp(1,1) * pha + lp(1,2);
% 
% elseif( pha>=pt(2) && pha<=pt(3) )
%     I = smooth_piecewise(pha+params.stim_start, params.t_amp_top_real, params.k0_Herm,0,amp,params.y1_Herm);
% 
% elseif( pha>pt(3) && pha<=pt(4) )
%     I = lp(2,1) * pha + lp(2,2);
%     I = smooth_piecewise(pha+params.stim_start, params.t_amp_top_real+2*params.t_Herm_dt,0,params.k1_Herm,params.y1_Herm,amp);
% elseif( pha>pt(4) && pha<=pt(5) )
%     I = lp(2,1) * pha + lp(2,2);
% else
%     I = 0;
% end


% I=interp1(params.min_step_t, params.Stim_wave, te, 'previous', 0);

end

%[appendix]{"version":"1.0"}
%---
