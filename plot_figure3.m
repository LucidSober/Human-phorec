% plot_figure3.m
function plot_figure3(light_states, dark_states, params)
% Reproduces Figure 3 from Muangkram et al., 2023.
% This figure shows ATP expenditure in dark vs. light steady-state conditions.

    figure('Name', 'Figure 3 Reproduction', 'Position', [100, 100, 1200, 600]);

    % --- Get Dark Steady-State Data ---
    dark_currents = calculate_currents(dark_states, params);
    [~, ~, ~, dark_atp_photo] = calculate_ATP_consumption_accurate(dark_currents, dark_states, params, 0);
    atp_dark_components = calculate_atp_components(dark_currents, params);
    
    % --- Get Light Steady-State Data ---
    % --- FIX: The complete light_states struct is now passed directly to this function ---
    light_currents = calculate_currents(light_states, params);
    stim_light = params.light_intensities(end) * 0.43;
    [~, ~, ~, light_atp_photo] = calculate_ATP_consumption_accurate(light_currents, light_states, params, stim_light);
    atp_light_components = calculate_atp_components(light_currents, params);

    % --- Assemble Data for Plotting ---
    labels = {'I_{CNG}', 'I_{NCKX}', 'I_{CaL}', 'I_h', 'I_{NCX}', 'J_{NKCC1}', 'I_{L,Na}', 'I_{PMCA}', 'Photo-transduction', 'Total'};
    dark_values = zeros(1, length(labels));
    light_values = zeros(1, length(labels));
    
    dark_values(1) = atp_dark_components.ICNG;
    dark_values(2) = atp_dark_components.INCKX;
    dark_values(3) = atp_dark_components.ICaL;
    dark_values(4) = atp_dark_components.Ih;
    dark_values(5) = atp_dark_components.INCX;
    dark_values(6) = atp_dark_components.JNKCC1;
    dark_values(7) = atp_dark_components.IL_Na;
    dark_values(8) = atp_dark_components.IPMCA;
    dark_values(9) = dark_atp_photo;
    dark_values(10) = sum(dark_values(1:9));
    
    light_values(1) = atp_light_components.ICNG;
    light_values(2) = atp_light_components.INCKX;
    light_values(3) = atp_light_components.ICaL;
    light_values(4) = atp_light_components.Ih;
    light_values(5) = atp_light_components.INCX;
    light_values(6) = atp_light_components.JNKCC1;
    light_values(7) = atp_light_components.IL_Na;
    light_values(8) = atp_light_components.IPMCA;
    light_values(9) = light_atp_photo;
    light_values(10) = sum(light_values(1:9));

    % --- Create the Bar Chart ---
    bar_data = [dark_values; light_values]';
    b = bar(bar_data, 'grouped');
    b(1).FaceColor = 'k'; % Dark bars are black
    b(2).FaceColor = 'w'; % Light bars are white
    
    set(gca, 'XTickLabel', labels, 'FontSize', 12);
    ylabel('ATP expenditure (molecule s^{-1})', 'FontSize', 14);
    ylim([0 8e7]);
    title('Figure 3: ATP Expenditure in Dark (D) vs. Light (L)', 'FontSize', 16);
    legend({'Dark', 'Light'});
    grid on;
    box on;
end