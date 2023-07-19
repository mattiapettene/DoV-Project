function simulation = getSimulationParams()
 
    % Timings
    times.t0        = 0;     % [s]  <--- starting time
    times.step_size = 1e-4; % [s]  <--- discrete solver step
    
    
    simulation.times = times;

end