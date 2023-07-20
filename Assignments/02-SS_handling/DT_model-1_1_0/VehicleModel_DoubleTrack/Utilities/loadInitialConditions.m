function ICs = loadInitialConditions(u0)

    % ----------------------------------------------------------------
    %% Function purpose: define the initial conditions for vehicle states
    % ----------------------------------------------------------------

    %% Load vehicle and auxiliary data
    vehicle_data = getVehicleDataStruct();
    m  = vehicle_data.vehicle.m; 
    Lr = vehicle_data.vehicle.Lr; 
    Lf = vehicle_data.vehicle.Lf;   
    L  = Lf + Lr;
    Rf = vehicle_data.front_wheel.Rf;
    Rr = vehicle_data.rear_wheel.Rr; 
    CAzf = vehicle_data.aerodynamics.CAzf;  
    CAzr = vehicle_data.aerodynamics.CAzr;
    g  = vehicle_data.vehicle.g; 
    
    
    %% ICs for all the states
    x0         = 0;           % [m] initial x-coordinate of the vehicle CoM
    y0         = 0;           % [m] initial y-coordinate of the vehicle CoM
    psi0       = 0;           % [rad] initial yaw angle
    %u0         = 100/3.6;      % [m/s] initial forward speed
    v0         = 0;           % [m/s] initial lateral speed
    Omega0     = 0;           % [rad/s] initial yaw rate
    Fz_rr0     = m*g*Lf/(2*L) + (CAzr*u0^2)/2;   % [N] initial vertical load for the rear right wheel 
    Fz_rl0     = m*g*Lf/(2*L) + (CAzr*u0^2)/2;   % [N] initial vertical load for the rear left wheel 
    Fz_fr0     = m*g*Lr/(2*L) + (CAzf*u0^2)/2;   % [N] initial vertical load for the front right wheel 
    Fz_fl0     = m*g*Lr/(2*L) + (CAzf*u0^2)/2;   % [N] initial vertical load for the front left wheel 
    delta0     = 0;           % [rad] initial steering angle
    omega__rr0 = u0/Rr;       % [rad/s] initial rotational speed of rear right wheel
    omega__rl0 = u0/Rr;       % [rad/s] initial rotational speed of rear left wheel
    omega__fr0 = u0/Rf;       % [rad/s] initial rotational speed of front right wheel
    omega__fl0 = u0/Rf;       % [rad/s] initial rotational speed of front left wheel
    alpha__rr0 = 0;           % [rad] initial slip angle for rear right wheel
    alpha__rl0 = 0;           % [rad] initial slip angle for rear left wheel
    alpha__fr0 = 0;           % [rad] initial slip angle for front right wheel
    alpha__fl0 = 0;           % [rad] initial slip angle for front left wheel
    kappa__rr0 = 0;           % [-] initial rear right wheel longitudinal slip
    kappa__rl0 = 0;           % [-] initial rear left wheel longitudinal slip
    kappa__fr0 = 0;           % [-] initial front right wheel longitudinal slip
    kappa__fl0 = 0;           % [-] initial front left wheel longitudinal slip
    
    ICs = [x0, y0, psi0, u0, v0, Omega0, ...
           Fz_rr0, Fz_rl0, Fz_fr0, Fz_fl0, delta0, ...
           omega__rr0, omega__rl0, omega__fr0, omega__fl0, ...
           alpha__rr0, alpha__rl0, alpha__fr0, alpha__fl0, ...
           kappa__rr0, kappa__rl0, kappa__fr0, kappa__fl0];
       
end

