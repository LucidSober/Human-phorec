% plot_figure5.m
function plot_figure5(results_constant_light, dark_states, light_intensities)
% Reproduces Figure 5 from Muangkram et al., 2023.
% Shows time courses and steady-state values vs. light intensity for constant light simulations.

    figure('Name', 'Figure 5 Reproduction', 'Position', [150, 150, 1000, 700]);
    colors = gray(length(light_intensities) + 2); % Black to gray colormap
    colors = colors(1:end-2,:);

    % --- Subplot A: Membrane Potential vs. Time ---
    subplot(2, 3, 1); hold on;
    for i = 1:length(light_intensities)
        plot(results_constant_light{i}.t - 1, results_constant_light{i}.Vm, 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    title('(A) Membrane Potential');
    xlabel('Time (s)');
    ylabel('(mV)');
    xlim([-1 8]);
    grid on; box on;

    % --- Subplot B: Total ATP Expenditure vs. Time ---
    subplot(2, 3, 2); hold on;
    for i = 1:length(light_intensities)
        plot(results_constant_light{i}.t - 1, results_constant_light{i}.ATP_total, 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    title('(B) Total Energy Expenditure');
    xlabel('Time (s)');
    ylabel('(molecules s^{-1})');
    xlim([-1 8]);
    grid on; box on;

    % --- Extract endpoint data for plots C, D, E, F ---
    end_values = struct();
    end_time = 5; % As specified in paper
    
    num_intensities = length(light_intensities);
    end_values.ATP_photo = zeros(num_intensities, 1);
    end_values.ATP_total = zeros(num_intensities, 1);
    end_values.ATP_Na = zeros(num_intensities, 1);
    end_values.ATP_Ca = zeros(num_intensities, 1);
    end_values.K = zeros(num_intensities, 1);
    end_values.Cl = zeros(num_intensities, 1);
    end_values.Na = zeros(num_intensities, 1);
    end_values.Caos = zeros(num_intensities, 1);
    end_values.Cais = zeros(num_intensities, 1);

    for i = 1:num_intensities
        res = results_constant_light{i};
        idx = find(res.t >= end_time, 1);
        if isempty(idx), idx = length(res.t); end
        
        end_values.ATP_photo(i) = res.ATP_cascade(idx);
        end_values.ATP_total(i) = res.ATP_total(idx);
        end_values.ATP_Na(i) = res.ATP_Na(idx);
        end_values.ATP_Ca(i) = res.ATP_Ca(idx);
        end_values.K(i) = res.Ki(idx);
        end_values.Cl(i) = res.Cli(idx);
        end_values.Na(i) = res.Nai(idx);
        end_values.Caos(i) = res.Caos(idx);
        end_values.Cais(i) = res.Cais(idx);
    end
    
    % --- Subplot C: Phototransduction ATP vs. Light Intensity ---
    subplot(2, 3, 3);
    semilogx(light_intensities, end_values.ATP_photo, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1);
    title('(C) Phototransduction Energy');
    xlabel('J_{hv} (photons \mu m^{-2} s^{-1})');
    ylabel('(molecules s^{-1})');
    grid on; box on;
    
    % --- Subplot D: ATP Breakdown vs. Light Intensity ---
    subplot(2, 3, 4); hold on;
    % --- MODIFIED: Using plot + set(gca, 'XScale', 'log') for robust logarithmic scaling ---
    plot(light_intensities, end_values.ATP_total, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(light_intensities, end_values.ATP_Na, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    plot(light_intensities, end_values.ATP_Ca, 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1);
    set(gca, 'XScale', 'log');
    title('(D) Energy Expenditure Breakdown');
    xlabel('J_{hv} (photons \mu m^{-2} s^{-1})');
    ylabel('(molecules s^{-1})');
    legend({'Total', 'Required for Na^+', 'Required for Ca^{2+}'}, 'Location', 'southwest');
    ylim([0 8e7]);
    xlim([1 1e5]);
    set(gca, 'XTick', [1e0, 1e1, 1e2, 1e3, 1e4, 1e5]);
    grid on; box on;

    % --- Subplot E: K+ and Cl- Alterations vs. Light Intensity ---
    alt_K = 100 * (end_values.K - dark_states.Ki) / dark_states.Ki;
    alt_Cl = 100 * (end_values.Cl - dark_states.Cli) / dark_states.Cli;
    subplot(2, 3, 5); hold on;
    % --- MODIFIED: Using plot + set(gca, 'XScale', 'log') for robust logarithmic scaling ---
    plot(light_intensities, alt_K, 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(light_intensities, alt_Cl, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
    set(gca, 'XScale', 'log');
    title('(E) Alterations of [K^+] and [Cl^-]');
    xlabel('J_{hv} (photons \mu m^{-2} s^{-1})');
    ylabel('Alteration from initial (%)');
    legend({'[K^+]', '[Cl^-]'}, 'Location', 'southeast');
    ylim([-1.5 3.0]);
    xlim([1 1e5]);
    set(gca, 'XTick', [1e0, 1e1, 1e2, 1e3, 1e4, 1e5]);
    grid on; box on;
    
    % --- Subplot F: Na+ and Ca2+ Alterations vs. Light Intensity ---
    alt_Na = 100 * (end_values.Na - dark_states.Nai) / dark_states.Nai;
    alt_Caos = 100 * (end_values.Caos - dark_states.Caos) / dark_states.Caos;
    alt_Cais = 100 * (end_values.Cais - dark_states.Cais) / dark_states.Cais;
    subplot(2, 3, 6); hold on;
    % --- MODIFIED: Using plot + set(gca, 'XScale', 'log') for robust logarithmic scaling ---
    plot(light_intensities, alt_Na, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(light_intensities, alt_Caos, 'kd', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    plot(light_intensities, alt_Cais, 'x', 'MarkerSize', 10, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
    set(gca, 'XScale', 'log');
    title('(F) Alterations of [Na^+] and [Ca^{2+}]');
    xlabel('J_{hv} (photons \mu m^{-2} s^{-1})');
    ylabel('Alteration from initial (%)');
    legend({'[Na^+]', '[Ca^{2+}]_{os}', '[Ca^{2+}]_{is}'}, 'Location', 'southwest');
    ylim([-100 0]);
    xlim([1 1e5]);
    set(gca, 'XTick', [1e0, 1e1, 1e2, 1e3, 1e4, 1e5]);
    grid on; box on;

end