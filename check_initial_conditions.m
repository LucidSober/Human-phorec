% check_initial_conditions.m
function check_initial_conditions()
    % Checks the reasonableness of the initial conditions.
    params = init_parameters();
    states = init_state_variables();
    
    % Calculate initial currents
    currents = calculate_currents(states, params);
    
    fprintf('\nInitial Condition Check:\n');
    fprintf('Membrane Potential: %.2f mV\n', states.Vm);
    fprintf('Ion Concentrations:\n');
    fprintf('  Ki = %.1f mM, Nai = %.1f mM\n', states.Ki, states.Nai);
    fprintf('  Caos = %.3f uM, Cais = %.3f uM\n', states.Caos, states.Cais);
    fprintf('  cGMP = %.2f uM\n', states.cGMP);
    
    fprintf('\nInitial Currents:\n');
    fprintf('  ICNG = %.2f pA\n', currents.I_CNG);
    fprintf('  INCKX = %.2f pA\n', currents.I_NCKX);
    fprintf('  Total Dark Current (CNG+NCKX) = %.2f pA\n', currents.I_CNG + currents.I_NCKX);
    
    % Check the sum of phototransduction components
    R_total_check = states.R + states.R0 + states.R1 + states.R2 + states.R3 + ...
                    states.R4 + states.R5 + states.R6 + states.R_Gt;
    fprintf('\nPhototransduction Check:\n');
    fprintf('  Sum of R species = %.2e (should be close to Rtot = %.2e)\n', R_total_check, params.Rtot);
    
    if abs(currents.I_CNG + currents.I_NCKX + 15) > 5
        warning('Initial dark current deviates significantly from the expected value of ~-15 pA!');
    end
end