function simulation = getSimulationParams()
 
    % Timings
    times.t0        = 0;     % [s]  <--- starting time
    times.step_size = 0.5e-4; % [s]  <--- discrete solver step
    times.tf        = 200;    % [s]  <--- stop simulation time
    
    simulation.times = times;

end