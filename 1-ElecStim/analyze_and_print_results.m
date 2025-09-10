%[text] 此函数的简短摘要。
%[text] 此函数的详细说明。
function analyze_and_print_results(results, t_results, light_intensities)
    % 定义要输出的时间点和变量
    global params

    time_points = [0.5, 0.52, 0.7, 2.0];
    Cm = params.Cm;

    % 定义要打印的24个变量及其单位和转换因子
    vars_to_print = {
        'Vm', 'mV', 1;
        'ICNG', 'pA/pF', 1/Cm; 'INCKX', 'pA/pF', 1/Cm; 'INaK', 'pA/pF', 1/Cm;
        'Ih', 'pA/pF', 1/Cm; 'IKv', 'pA/pF', 1/Cm; 'IKCa', 'pA/pF', 1/Cm;
        'IClCa', 'pA/pF', 1/Cm; 'ICaL', 'pA/pF', 1/Cm; 'IPMCA', 'pA/pF', 1/Cm;
        'INCX', 'fA/pF', 1000/Cm; 'JNKCC1', 'mM/s', 1; 'JKCC2', 'mM/s', 1;
        'IL_K', 'fA/pF', 1000/Cm; 'IL_Na', 'fA/pF', 1000/Cm; 'IL_Cl', 'fA/pF', 1000/Cm;
        'IL_Caos', 'fA/pF', 1000/Cm; 'IL_Cais', 'fA/pF', 1000/Cm;
        'Caos', 'uM', 1; 'Cais', 'uM', 1; 'cGMP', 'uM', 1;
        'Ki', 'mM', 1; 'Nai', 'mM', 1; 'Cli', 'mM', 1;
    };

    for i = 1:length(results)
        fprintf('\n======================================================\n');
        fprintf('  结果分析: 光照强度 %.1f photons/μm²/s\n', light_intensities(i));
        fprintf('======================================================\n');
        
        for t_target = time_points
            % 找到最接近目标时间点的数据索引
            [~, idx] = min(abs(t_results{i} - t_target));
            actual_time = t_results{i}(idx);
            
            fprintf('\n  @ 时间点: %.4f s\n', actual_time);
            fprintf('  ------------------------------------\n');
            
            for v = 1:size(vars_to_print, 1)
                var_name = vars_to_print{v, 1};
                var_unit = vars_to_print{v, 2};
                var_conv = vars_to_print{v, 3};
                
                % 检查results中是否存在该字段
                if isfield(results{i}, var_name)
                    value = results{i}.(var_name)(idx) * var_conv;
                    fprintf('    - %-12s : %12.4f %s\n', var_name, value, var_unit);
                end
            end
        end
    end



    
end

%[appendix]{"version":"1.0"}
%---
