% ---------------------------------------------------------------
%% Main script for a basic simulation framework with a double track vehcile model
%  authors: 
%  rev. 1.0 Mattia Piccinini & Gastone Pietro Papini Rosati
%  rev. 2.0 Edoardo Pagot
%  date:
%  rev 1.0:    13/10/2020
%  rev 2.0:    16/05/2022
%  rev 2.1:    08/07/2022 (Biral)
%       - added Fz saturation. Correceted error in Fx
%       - initial condition is now parametric in initial speed
%       - changed the braking torque parameters to adapt to a GP2 model
% ----------------------------------------------------------------

% ----------------------------
%% Initialization
% ----------------------------
initialize_environment;

% ----------------------------
%% Load vehicle data
% ----------------------------

% test_tyre_model; % some plot to visualize the curvers resulting from the
% loaded data

vehicle_data = getVehicleDataStruct();
% pacejkaParam = loadPacejkaParam();

% ----------------------------
%% Define initial conditions for the simulation and type of test
% ----------------------------
% Description:
% - the simulation time is increased: Tf = 150s in order to see better results.
% - if switch_test_type = 1 -> speed ramp test at constant steering angle
% - if switch_test_type = 2 -> steer ramp test at constant forward speed
% Parameters description for each test:
% - t1 = time at which the transient ends, during t1 -> the desired
%           velocity is reached (thanks to PID);
% - t2 = time at which the ramp (during ramp test) ends -> after t2 the
%           imposed velocity/steer remains constant up to Tf
% - ax_imposed = longitudinal acceleration (low value) imposed during the
%           speed ramp test
% - steer_angle_slope = the gradient of the steer angle during steer ramp test is imposed equal
%           to 0.2 (in order to achive 20 deg of steer angle at the end of
%           simulation)

V0 = 50/3.6; % Initial speed
X0 = loadInitialConditions(V0);

V_init = V0;
ax_imposed = 0.12; %m/s^2
Tf = 150;
t1 = Tf/10;
t2 = Tf - Tf/10;
steer_angle_slope = 0.2;

const_steer_angle = 10; % [deg]
const_v_des = 50/3.6; % [m/s]

switch_test_type = 1; %1 = speed ramp test with const steer, 2 = steer ramp test  with const speed;

% ----------------------------
%% Simulation parameters
% ----------------------------
simulationPars = getSimulationParams(); 
Ts = simulationPars.times.step_size;  % integration step for the simulation (fixed step)
T0 = simulationPars.times.t0;         % starting time of the simulation
Tf = simulationPars.times.tf;         % stop time of the simulation

% ----------------------------
%% Start Simulation
% ----------------------------
fprintf('Starting Simulation\n')
tic;
model_sim = sim('Vehicle_Model_2Track');
elapsed_time_simulation = toc;
fprintf('Simulation completed\n')
fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

% ----------------------------
%% Post-Processing
% ----------------------------
dataAnalysis(model_sim,vehicle_data,Ts);
vehicleAnimation(model_sim,vehicle_data,Ts);
