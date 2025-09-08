% calculate_atp_components.m
function atp_breakdown = calculate_atp_components(currents, params)
% Calculates the ATP expenditure for each component causing Na+ or Ca2+ influx.
% This is based on the methodology described in the paper (Muangkram et al., 2023).

    NA = params.NA;      % Avogadro's constant
    F = params.F;        % Faraday's constant
    V_is = params.V_is;  % Inner segment volume
    
    atp_breakdown = struct();

    % --- ATP for Na+ Extrusion, broken down by source of Na+ influx ---
    % 1 ATP is used to pump out 3 Na+ ions.
    
    % Convert pA to A (C/s)
    I_CNG_Na_A = currents.I_CNG_Na * 1e-12;
    I_h_Na_A = currents.I_h_Na * 1e-12;
    I_CaL_Na_A = currents.I_CaL_Na * 1e-12;
    I_L_Na_A = currents.I_L_Na * 1e-12;
    
    % Stoichiometry of Na+ ions per cycle for exchangers
    I_NCKX_A = currents.I_NCKX * 1e-12; % 4 Na+ in per cycle
    I_NCX_A = currents.I_NCX * 1e-12;   % 3 Na+ in per cycle
    
    % Note: Current is defined as positive charge flowing out. Influx is negative current.
    % Moles of Na+ per second = -Current / (z * F), where z=1 for Na+
    % Molecules per second = Moles/sec * NA
    % ATP molecules per second = (Na+ molecules/sec) / 3
    
    atp_breakdown.ICNG = (-I_CNG_Na_A / F) * NA / 3;
    atp_breakdown.Ih = (-I_h_Na_A / F) * NA / 3;
    atp_breakdown.ICaL = (-I_CaL_Na_A / F) * NA / 3;
    atp_breakdown.IL_Na = (-I_L_Na_A / F) * NA / 3;
    
    % For exchangers, the net current reflects charge movement, not just Na+ movement.
    % NCKX: 4 Na+ in, 1 Ca2+ out, 1 K+ out. Net charge in = 4 - 2 - 1 = +1.
    % So, cycles/sec = I_NCKX_A / (1 * 1.602e-19)
    % Na+ ions/sec = 4 * cycles/sec
    charge_e = 1.60217663e-19;
    nckx_cycles_per_sec = -I_NCKX_A / charge_e; % Negative current is influx of net charge
    atp_breakdown.INCKX = (4 * nckx_cycles_per_sec) / 3;

    % NCX: 3 Na+ in, 1 Ca2+ out. Net charge in = 3 - 2 = +1.
    ncx_cycles_per_sec = -I_NCX_A / charge_e;
    atp_breakdown.INCX = (3 * ncx_cycles_per_sec) / 3;
    
    % For J_NKCC1 (electroneutral), we assume it's negligible for ATP cost as in the paper's Fig 3.
    atp_breakdown.JNKCC1 = 0;

    % --- ATP for Ca2+ Extrusion (via IPMCA) ---
    % 1 ATP is used to pump out 1 Ca2+ ion.
    % The paper attributes the I_PMCA cost directly.
    I_PMCA_A = currents.I_PMCA * 1e-12;
    % Moles of Ca2+ per second = Current / (z * F), where z=2
    % Molecules per second = Moles/sec * NA
    % ATP molecules per second = Ca2+ molecules/sec / 1
    atp_breakdown.IPMCA = (I_PMCA_A / (2 * F)) * NA;
    
    % Set any negative (computationally noisy) values to zero
    fields = fieldnames(atp_breakdown);
    for i = 1:length(fields)
        if atp_breakdown.(fields{i}) < 0
            atp_breakdown.(fields{i}) = 0;
        end
    end
end