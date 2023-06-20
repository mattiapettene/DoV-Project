function simulation = getSimulationParams()
 
    % Timings
    times.t0        = 0;     % [s]  <--- starting time
    times.step_size = 1e-4; % [s]  <--- discrete solver step
    times.tf        = 150;    % [s]  <--- stop simulation time
    
    simulation.times = times;

end