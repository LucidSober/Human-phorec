%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function plot_final_results_figure(results, t_results,t_range)
% PLOT_FINAL_RESULTS_FIGURE2 最终版的绘图函数
%   - 完整复现论文中的 Figure 2
%   - 使用彩色曲线以便于区分
    global params

    figure('Name', 'Figure 2 Reproduction (Color)', 'Position', [50, 50, 1400, 1000]);
    
    % %%%%%%%%%%%%%%%%%%%%%% 关键修改 %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 将颜色方案从灰度(gray)修改为彩色(jet)
    % colors = jet(length(results));
    colors = [linspace(0,1,length(results))' ...   % R 通道：从 0 到 1
          zeros(length(results),1) ...         % G 通道：保持 0
          linspace(1,0,length(results))'];     % B 通道：从 1 到 0
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Cm = params.Cm; 
    time_range = t_range;

    % --- 辅助绘图函数 ---
    function plot_subplot(pos, data_cell, ylabel_txt, title_txt, y_limits, unit_conv)
        subplot(6, 4, pos);
        hold on;
        for i = 1:length(results)
            plot(t_results{i}, data_cell{i} * unit_conv, 'Color', colors(i,:), 'LineWidth', 1.2);
        end
        title(title_txt, 'FontSize', 10);
        ylabel(ylabel_txt);
        xlim(time_range);
        ylim(y_limits);
        grid on;
        box on;
        set(gca, 'FontSize', 8);
    end

    %% 第1行
    plot_subplot(1, cellfun(@(x) x.Vm, results, 'UniformOutput', false), '(mV)', '1.Membrane potential', [-60, -30], 1);
    plot_subplot(2, cellfun(@(x) x.ICNG, results, 'UniformOutput', false), '(pA/pF)', '2.I_{CNG}_OS', [-5.0, 1.0], 1/Cm);
    plot_subplot(3, cellfun(@(x) x.INCKX, results, 'UniformOutput', false), '(pA/pF)', '3.I_{NCKX_OS}', [-0.30, 0.10], 1/Cm);
    plot_subplot(4, cellfun(@(x) x.INaK, results, 'UniformOutput', false), '(pA/pF)', '4.I_{NaK}', [1.80, 2.40], 1/Cm);

    %% 第2行
    plot_subplot(5, cellfun(@(x) x.Ih, results, 'UniformOutput', false), '(pA/pF)', '5.I_h', [-5.0, 1.0], 1/Cm);
    plot_subplot(6, cellfun(@(x) x.IKv, results, 'UniformOutput', false), '(pA/pF)', '6.I_{Kv}', [0.00, 3.0], 1/Cm);
    plot_subplot(7, cellfun(@(x) x.IKCa, results, 'UniformOutput', false), '(pA/pF)', '7.I_{KCa}', [-0.02, 0.10], 1/Cm);
    plot_subplot(8, cellfun(@(x) x.IClCa, results, 'UniformOutput', false), '(pA/pF)', '8.I_{ClCa}', [-0.50, 0.10], 1/Cm);

    %% 第3行
    plot_subplot(9, cellfun(@(x) x.ICaL, results, 'UniformOutput', false), '(pA/pF)', '9.I_{CaL}', [-2.0, 0.5], 1/Cm);
    plot_subplot(10, cellfun(@(x) x.IPMCA, results, 'UniformOutput', false), '(pA/pF)', '10.I_{PMCA}', [0.0, 2.5], 1/Cm);
    plot_subplot(11, cellfun(@(x) x.INCX, results, 'UniformOutput', false), '(fA/pF)', '11.I_{NCX}', [-15, 5], 1000/Cm);
    % 绘图 J_NKCC1，单位从 mM 转换为 µM
    plot_subplot(12, cellfun(@(x) x.JNKCC1 * -1000, results, 'UniformOutput', false), '(\muM)', '12.J_{NKCC1}', [23, 25], 1);
%% 第4行
% 绘图 J_KCC2，单位从 mM 转换为 µM
    plot_subplot(13, cellfun(@(x) x.JKCC2 * 1000, results, 'UniformOutput', false), '(\muM)', '13.J_{KCC2}', [-28.14, -28.015], 1);
    plot_subplot(14, cellfun(@(x) x.IL_K, results, 'UniformOutput', false), '(fA/pF)', '14.I_{L,K}', [0.80, 2.00], 1000/Cm);
    plot_subplot(15, cellfun(@(x) x.IL_Na, results, 'UniformOutput', false), '(fA/pF)', '15.I_{L,Na}', [-3.10, -2.20], 1000/Cm);
    plot_subplot(16, cellfun(@(x) x.IL_Cl, results, 'UniformOutput', false), '(fA/pF)', '16.I_{L,Cl}', [-3.0, 1.0], 1000/Cm);

    %% 第5行
    plot_subplot(17, cellfun(@(x) x.IL_Caos, results, 'UniformOutput', false), '(fA/pF)', '17.I_{L,Caos_OS}', [-1.20, -0.80], 1000/Cm);
    plot_subplot(18, cellfun(@(x) x.IL_Cais, results, 'UniformOutput', false), '(fA/pF)', '18.I_{L,Cais}', [-14.5, -10.5], 1000/Cm);
    plot_subplot(19, cellfun(@(x) x.Caos, results, 'UniformOutput', false), '(\muM)', '[19.Ca^{2+}]_{os}', [0.000, 0.120], 1);
    plot_subplot(20, cellfun(@(x) x.Cais, results, 'UniformOutput', false), '(\muM)', '[20.Ca^{2+}]_{is}', [0.024, 0.072], 1);
    
    %% 第6行
    plot_subplot(21, cellfun(@(x) x.cGMP, results, 'UniformOutput', false), '(\muM)', '[21.cGMP]', [0.0, 8.0], 1);
    plot_subplot(22, cellfun(@(x) x.Ki, results, 'UniformOutput', false), '(mM)', '[22.K^+]_i', [139.8, 141.3], 1);
    plot_subplot(23, cellfun(@(x) x.Nai, results, 'UniformOutput', false), '(mM)', '[23.Na^+]_i', [3.4, 5.4], 1);
    plot_subplot(24, cellfun(@(x) x.Cli, results, 'UniformOutput', false), '(mM)', '[24.Cl^-]_i', [29.84, 30.04], 1);

    % 为最后一行添加X轴标签
    for k = 21:24
        xlabel(subplot(6, 4, k), 'Time (s)');
    end
    
    sgtitle('Figure 2 Reproduction: Electrophysiological Responses', 'FontSize', 16, 'FontWeight', 'bold');

    ax = findall(gcf, 'Type', 'axes');        % 获取当前图中所有 axes 对象
    set(ax, 'YLimMode', 'auto'); 
    set(ax, 'YLimitMethod', 'padded');
end

%[appendix]{"version":"1.0"}
%---
