

### **Rod Photoreceptor Model Run Instructions**

This guide will assist you in successfully running and using the provided MATLAB rod photoreceptor model.

#### 1. Prerequisites

Before running the model, please ensure you have completed the following preparations:

  File Preparation: Ensure that all the following `.m` files are located in the same folder. This is essential for the model to run correctly:

      main.m
      derivatives.m
      init_parameters.m
      init_state_variables.m
      calculate_currents.m
      phototransduction_cascade.m
      plot_final_results_figure2.m
      calculate_ATP_consumption_accurate.m

    Helper Script (Optional):

       check_initial_conditions.m`: Can be used to verify the reasonability of the initial parameter settings.
      
#### 2. Running the Simulation

1.  Open MATLAB.
2.  Navigate to the File Directory: In the MATLAB interface, use the address bar or the `cd` command in the Command Window to change the current working directory to the folder containing all the `.m` files.
    ```matlab
    % Example: cd 'C:\Your\Folder\Path\To\Model_Files'
    cd 'path\to\your\model_files' 
    ```
3.  Execute the Main Script: In the MATLAB Command Window, type the name of the main script, `main`, and press the **Enter** key.
    or open the main script manually and click the green button to run.
	

The script will begin execution automatically. No further actions are required.

#### 3. Expected Output

As the script runs, you will observe the following outputs:

1.  **Command Window Output**:

       The script will first print progress messages about initialization and finding the steady state.
        ```
        1. Initializing parameters and state variables...
        2. Finding the steady state in dark conditions...
        ```
      This will be followed by a "Steady-State Check," which compares the model's calculated dark-adapted steady-state values against expected values.
      The script will then begin simulating the response to different light intensities and will print its progress.
        ```
        3. Starting light response simulations...
           Simulating 1/8: Light Intensity = 1.7 photons/um^2/s
           ...
           Simulating 8/8: Light Intensity = 4630.0 photons/um^2/s
        ```
      After the simulations are complete, a detailed analysis of key variables at specific time points (e.g., 0.5s, 0.52s) will be printed.
      Finally, messages indicating that plotting and the simulation are complete will appear.
        ```
        5. Plotting results...
        Simulation complete!
        ```

2.  Graphical Output:

       Upon completion, a MATLAB figure window will automatically open with the title "Figure 2 Reproduction: Electrophysiological Responses to Light Flashes".
      This figure contains a **6x4 grid of 24 subplots**, showing the time courses of various electrophysiological and ionic concentration variables (e.g., membrane potential, various ion currents, calcium concentrations) for each simulated light intensity.

#### 4. Customization Guide

You can easily modify key parameters to suit your research needs. All recommended modifications should be made in the `main.m` file.

  Modifying Light Intensity：

      In `main.m`, locate the following variable:
        ```matlab
        light_intensities = [1.7, 4.8, 15.2, 39.4, 125, 444, 1406, 4630];
        ```
      You can modify the values in this array or add/remove elements to simulate the light intensities you are interested in (units: photons/μm²/s).

  Modifying Simulation Time & Flash Parameters:

      In `main.m`, locate the following variables:
        ```matlab
        t_total = 5.0; % s
        flash_start = 0.5; % s
        flash_duration = 0.02; % s
        ```
      `t_total`: The total simulation duration.
      `flash_start`: The onset time of the light flash.
      `flash_duration`: The duration of the light flash.

  Modifying Plotting Ranges:

      Open the `plot_final_results_figure2.m` file.
      The `plot_subplot` function calls within this file are responsible for drawing each subplot. Its fifth argument, `y_limits`, defines the y-axis range.
      For example, to modify the plotting range for `IKCa`, find the corresponding line:
        ```matlab
        % Original line:
        % plot_subplot(7, cellfun(@(x) x.IKCa, results, 'UniformOutput', false), '(pA/pF)', 'I_{KCa}', [-0.02, 0.22], 1/Cm);

        % To change the range to, for example, [-0.05, 0.25], modify it as follows:
        plot_subplot(7, cellfun(@(x) x.IKCa, results, 'UniformOutput', false), '(pA/pF)', 'I_{KCa}', [-0.05, 0.25], 1/Cm);
        ```

#### **5. File Descriptions**

  main.m`: **Main Script**. Sets up parameters, runs the simulation, and calls the plotting function.
  derivatives.m`: **Core ODE Function**. Defines the differential equations for all 94 state variables, which are called by the `ode15s` solver.
  init_parameters.m`: **Parameter Initialization**. Defines all physical and chemical constants used in the model.
  init_state_variables.m`: **State Variable Initialization**. Sets the dark-adapted initial values for all state variables.
  calculate_currents.m`: **Current Calculation Function**. Calculates all transmembrane ionic currents based on the current state of the model.
  phototransduction_cascade.m`: **Phototransduction Cascade Function**. Simulates the biochemical reactions from photon absorption to PDE activation.
  plot_final_results_figure2.m`: **Plotting Function**. Visualizes the simulation results.
  calculate_ATP_consumption_accurate.m`: **ATP Consumption Function**. Used for post-processing analysis of energy consumption.