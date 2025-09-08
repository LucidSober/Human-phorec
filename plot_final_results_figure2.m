% plot_final_results_figure2.m
function plot_final_results_figure2(results, t_results, params)
% PLOT_FINAL_RESULTS_FIGURE2 The final plotting function.
%   - Completely reproduces Figure 2 from the source paper.
%   - Uses a color map for better curve discrimination.

    figure('Name', 'Figure 2 Reproduction (Color)', 'Position', [50, 50, 1400, 1000]);
    
    colors = jet(length(results)); % Use 'jet' colormap for colored lines

    Cm = params.Cm; 
    time_range = [0, 5]; % x-axis time range in seconds

    % --- Plotting helper function ---
    function plot_subplot(pos, data_cell, ylabel_txt, title_txt, y_limits, unit_conv)
        subplot(6, 4, pos);
        hold on;
        for i = 1:length(results)
            plot(t_results{i}, data_cell{i} * unit_conv, 'Color', colors(i,:), 'LineWidth', 1.2);
        end
        title(title_txt, 'FontSize', 10);
        ylabel(ylabel_txt);
        xlim(time_range);
        if ~isempty(y_limits)
            ylim(y_limits);
        end
        grid on;
        box on;
        set(gca, 'FontSize', 8);
    end

    %% Row 1
    plot_subplot(1, cellfun(@(x) x.Vm, results, 'UniformOutput', false), '(mV)', 'Membrane Potential', [-60, -30], 1);
    plot_subplot(2, cellfun(@(x) x.ICNG, results, 'UniformOutput', false), '(pA/pF)', 'I_{CNG}', [-5.0, 1.0], 1/Cm);
    plot_subplot(3, cellfun(@(x) x.INCKX, results, 'UniformOutput', false), '(pA/pF)', 'I_{NCKX}', [-0.30, 0.10], 1/Cm);
    plot_subplot(4, cellfun(@(x) x.INaK, results, 'UniformOutput', false), '(pA/pF)', 'I_{NaK}', [1.80, 2.40], 1/Cm);

    %% Row 2
    plot_subplot(5, cellfun(@(x) x.Ih, results, 'UniformOutput', false), '(pA/pF)', 'I_h', [-5.0, 1.0], 1/Cm);
    plot_subplot(6, cellfun(@(x) x.IKv, results, 'UniformOutput', false), '(pA/pF)', 'I_{Kv}', [0.00, 3.0], 1/Cm);
    plot_subplot(7, cellfun(@(x) x.IKCa, results, 'UniformOutput', false), '(pA/pF)', 'I_{KCa}', [-0.02, 0.26], 1/Cm);
    
    plot_subplot(8, cellfun(@(x) x.IClCa, results, 'UniformOutput', false), '(pA/pF)', 'I_{ClCa}', [-0.50, 0.10], 1/Cm);

    %% Row 3
    plot_subplot(9, cellfun(@(x) x.ICaL, results, 'UniformOutput', false), '(pA/pF)', 'I_{CaL}', [-2.0, 0.5], 1/Cm);
    plot_subplot(10, cellfun(@(x) x.IPMCA, results, 'UniformOutput', false), '(pA/pF)', 'I_{PMCA}', [0.0, 2.5], 1/Cm);
    plot_subplot(11, cellfun(@(x) x.INCX, results, 'UniformOutput', false), '(fA/pF)', 'I_{NCX}', [-15, 5], 1000/Cm);
    plot_subplot(12, cellfun(@(x) x.JNKCC1 * 1000, results, 'UniformOutput', false), '(uM/s)', 'J_{NKCC1}', [23, 25], 1);

    %% Row 4
    plot_subplot(13, cellfun(@(x) x.JKCC2 * 1000, results, 'UniformOutput', false), '(uM/s)', 'J_{KCC2}', [-28.14, -28.015], 1);
    plot_subplot(14, cellfun(@(x) x.IL_K, results, 'UniformOutput', false), '(fA/pF)', 'I_{L,K}', [0.80, 2.00], 1000/Cm);
    plot_subplot(15, cellfun(@(x) x.IL_Na, results, 'UniformOutput', false), '(fA/pF)', 'I_{L,Na}', [-3.10, -2.20], 1000/Cm);
    plot_subplot(16, cellfun(@(x) x.IL_Cl, results, 'UniformOutput', false), '(fA/pF)', 'I_{L,Cl}', [-3.0, 1.0], 1000/Cm);

    %% Row 5
    plot_subplot(17, cellfun(@(x) x.IL_Caos, results, 'UniformOutput', false), '(fA/pF)', 'I_{L,Caos}', [-1.20, -0.80], 1000/Cm);
    plot_subplot(18, cellfun(@(x) x.IL_Cais, results, 'UniformOutput', false), '(fA/pF)', 'I_{L,Cais}', [-14.5, -10.5], 1000/Cm);
   
    plot_subplot(19, cellfun(@(x) x.Caos, results, 'UniformOutput', false), '(uM)', '[Ca^{2+}]_{os}', [0.000, 0.290], 1);
    
    plot_subplot(20, cellfun(@(x) x.Cais, results, 'UniformOutput', false), '(uM)', '[Ca^{2+}]_{is}', [0.024, 0.130], 1);
    
    %% Row 6
    plot_subplot(21, cellfun(@(x) x.cGMP, results, 'UniformOutput', false), '(uM)', '[cGMP]', [0.0, 8.0], 1);
    plot_subplot(22, cellfun(@(x) x.Ki, results, 'UniformOutput', false), '(mM)', '[K^+]_i', [139.8, 141.3], 1);
    plot_subplot(23, cellfun(@(x) x.Nai, results, 'UniformOutput', false), '(mM)', '[Na^+]_i', [3.4, 5.4], 1);
    plot_subplot(24, cellfun(@(x) x.Cli, results, 'UniformOutput', false), '(mM)', '[Cl^-]_i', [29.84, 30.04], 1);

    % Add X-axis labels to the bottom row of subplots
    for k = 21:24
        xlabel(subplot(6, 4, k), 'Time (s)');
    end
    
    sgtitle('Figure 2 Reproduction: Electrophysiological Responses to Light Flashes', 'FontSize', 16, 'FontWeight', 'bold');
end