function test_corrected_RC_stim()
% 测试修正的双段RC刺激函数并与原始cGMP数据对比
% 调用 corrected_RC_stim.m 函数

close all
clear all
clc

%% 1. 加载原始数据
try
    load cGMP
    cGMP = -cGMP;
    
    % 数据预处理
    amp = 40e-6;
    min_val = min(cGMP(:));
    max_val = max(cGMP(:));
    cGMP = (cGMP - min_val) / (max_val - min_val) * amp;
    cGMP_back = cGMP(976:end);
    
    % 时间轴
    t_start = 1.25;
    t_end = 5;
    t_original = linspace(t_start, t_end, length(cGMP_back))';
    
    fprintf('成功加载cGMP数据，数据点数: %d\n', length(cGMP_back));
    has_data = true;
    
catch
    fprintf('警告：未找到cGMP.mat文件\n');
    has_data = false;
end

%% 2. 生成拟合曲线
% 创建高分辨率时间轴用于拟合曲线
t_fit = linspace(1, 3, 1000)';

% 调用修正的RC刺激函数
I_fit = corrected_RC_stim(t_fit);

% 在原始数据时间点上计算拟合值（用于计算拟合优度）
if has_data
    I_predicted = corrected_RC_stim(t_original);
    
    % 计算拟合优度 (只在刺激时间段内)
    stim_mask = (t_original >= 1.25) & (t_original <= 2.5);
    t_stim = t_original(stim_mask);
    I_stim_original = cGMP_back(stim_mask);
    I_stim_predicted = I_predicted(stim_mask);
    
    % R²计算
    ss_res = sum((I_stim_original - I_stim_predicted).^2);
    ss_tot = sum((I_stim_original - mean(I_stim_original)).^2);
    r_squared = 1 - (ss_res / ss_tot);
    
    % RMSE计算
    rmse = sqrt(mean((I_stim_original - I_stim_predicted).^2));
    
    fprintf('\n=== 拟合质量评估 ===\n');
    fprintf('拟合优度 R² = %.4f\n', r_squared);
    fprintf('均方根误差 RMSE = %.2e A (%.2f μA)\n', rmse, rmse*1e6);
end

%% 3. 绘制对比图
figure('Position', [100, 100, 1200, 800]);

% 主对比图
subplot(2, 2, 1);
if has_data
    plot(t_original, cGMP_back*1e6, 'b-', 'LineWidth', 2, 'DisplayName', '原始cGMP数据');
    hold on;
end
plot(t_fit, I_fit*1e6, 'r-', 'LineWidth', 2, 'DisplayName', '修正双段RC拟合');

% 标记关键点
key_times = [1.25, 1.75, 2.25];
key_currents = corrected_RC_stim(key_times);
scatter(key_times, key_currents*1e6, 100, 'ro', 'filled', 'DisplayName', '关键点');

% 添加垂直线标记
xline(1.25, 'g--', 'LineWidth', 1, 'DisplayName', '起始点');
xline(1.75, 'k--', 'LineWidth', 1, 'DisplayName', '衔接点');
xline(2.25, 'm--', 'LineWidth', 1, 'DisplayName', '结束点');

xlabel('时间 (s)');
ylabel('电流 (μA)');
title('修正双段RC模型 vs 原始数据');
legend('Location', 'northeast');
grid on;
xlim([1.2, 2.8]);

if has_data
    % 添加拟合信息到图上
    text(1.3, 35, sprintf('R² = %.3f', r_squared), 'FontSize', 12, 'BackgroundColor', 'white');
end

% 残差分析（如果有原始数据）
if has_data
    subplot(2, 2, 2);
    residuals = (I_stim_original - I_stim_predicted) * 1e6;  % 转换为μA
    plot(t_stim, residuals, 'ko-', 'MarkerSize', 4, 'LineWidth', 1);
    xlabel('时间 (s)');
    ylabel('残差 (μA)');
    title('拟合残差分析');
    grid on;
    yline(0, 'r--', 'LineWidth', 1);
    
    % 显示残差统计
    mean_residual = mean(residuals);
    std_residual = std(residuals);
    text(min(t_stim) + 0.1, max(residuals)*0.8, ...
         sprintf('均值: %.2f μA\n标准差: %.2f μA', mean_residual, std_residual), ...
         'FontSize', 10, 'BackgroundColor', 'white');
end

% 分段显示
subplot(2, 2, 3);
t1 = linspace(1.25, 1.75, 100);
t2 = linspace(1.75, 2.25, 100);
I1 = corrected_RC_stim(t1);
I2 = corrected_RC_stim(t2);

plot(t1, I1*1e6, 'g-', 'LineWidth', 3, 'DisplayName', '第一段：偏移+充电');
hold on;
plot(t2, I2*1e6, 'b-', 'LineWidth', 3, 'DisplayName', '第二段：放电');

% 显示衔接点
plot(1.75, 20, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red', 'DisplayName', '衔接点');

xlabel('时间 (s)');
ylabel('电流 (μA)');
title('分段显示');
legend('Location', 'northeast');
grid on;

% 在衔接点添加文本说明
text(1.76, 22, '20 μA', 'FontSize', 10, 'BackgroundColor', 'yellow');

% 函数参数显示
subplot(2, 2, 4);
axis off;

% 获取函数参数（从函数文件中读取或者设定）
I_offset = 2e-5;
I_amp = 2e-5;
t_start = 1.25;
t_mid = 1.75;
t_end = 2.25;

% 显示参数表格
param_text = {
    '=== 函数参数 ==='
    ''
    sprintf('偏移量: %.0f μA', I_offset*1e6)
    sprintf('幅度: %.0f μA', I_amp*1e6)
    sprintf('起始时间: %.2f s', t_start)
    sprintf('衔接时间: %.2f s', t_mid)
    sprintf('结束时间: %.2f s', t_end)
    ''
    '=== 关键点 ==='
    ''
    sprintf('t=%.2fs: %.0f μA', t_start, corrected_RC_stim(t_start)*1e6)
    sprintf('t=%.2fs: %.0f μA', t_mid, corrected_RC_stim(t_mid)*1e6)
    sprintf('t=%.2fs: %.1f μA', t_end, corrected_RC_stim(t_end)*1e6)
    ''
    '=== 函数公式 ==='
    ''
    '第一段 (1.25s ≤ t ≤ 1.75s):'
    'I = 20 + 20×(1-exp(-(1.75-t)/τ₁))'
    ''
    '第二段 (t > 1.75s):'
    'I = 20×exp(-(t-1.75)/τ₂)'
};

text(0.05, 0.95, param_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontName', 'FixedWidth', 'FontSize', 10);

%% 4. 保存结果
if has_data
    % 保存对比数据
    comparison_data = [t_original, cGMP_back, I_predicted];
    writematrix(comparison_data, 'comparison_results.txt');
    fprintf('\n对比结果已保存到: comparison_results.txt\n');
end

% 保存拟合曲线数据
fit_data = [t_fit, I_fit];
writematrix(fit_data, 'fitted_curve_data.txt');
fprintf('拟合曲线数据已保存到: fitted_curve_data.txt\n');

%% 5. 输出总结信息
fprintf('\n=== 测试总结 ===\n');
fprintf('函数名: corrected_RC_stim.m\n');
fprintf('测试时间范围: %.2f - %.2f s\n', min(t_fit), max(t_fit));
fprintf('电流范围: %.1f - %.1f μA\n', min(I_fit)*1e6, max(I_fit)*1e6);
fprintf('衔接点: t=%.2fs, I=%.0fμA\n', 1.75, corrected_RC_stim(1.75)*1e6);

if has_data
    fprintf('数据拟合质量: R² = %.4f\n', r_squared);
    if r_squared > 0.95
        fprintf('拟合质量: 优秀\n');
    elseif r_squared > 0.90
        fprintf('拟合质量: 良好\n');
    elseif r_squared > 0.80
        fprintf('拟合质量: 一般\n');
    else
        fprintf('拟合质量: 需要改进\n');
    end
end

end
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
tau2 = 0.167;         % 第二段时间常数 (s)

% 初始化输出
I = zeros(size(t));

% 第一段：偏移 + 反向RC充电 (1.25s ≤ t ≤ 1.75s)
idx1 = (t >= t_start) & (t <= t_mid);
I(idx1) = I_offset + I_amp * (1 - exp(-(t_mid - t(idx1)) / tau1));

% 第二段：RC放电 (t > 1.75s)
idx2 = t > t_mid;
I(idx2) = I_offset * exp(-(t(idx2) - t_mid) / tau2);

end