%------------------------------------------------------------------------------
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
%  rev 2.2:    07/03/2023 (Taddei)
%  rev 3.0:    27/03/2023 (Rosati) Modified added the part of the planner
%  rev 3.2:    02/04/2023 (Rosati) Added the road_condition for each wheel
%                                  and clean variables and files
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%% Initialization
%------------------------------------------------------------------------------

initialize_environment;

%------------------------------------------------------------------------------
%% Load vehicle data
%------------------------------------------------------------------------------

% test_tyre_model; % some plot to visualize the curvers resulting from the
                   % loaded data

vehicle_data = getVehicleDataStruct();
% pacejkaParam = loadPacejkaParam();
L            = vehicle_data.vehicle.L;
tau_delta    = vehicle_data.steering_system.tau_D;
R0           = vehicle_data.tyre_data_f.R0;

% ----------------------------
%% Longitudinal controller parameters
% ----------------------------
longCTRparam = longitController();

% ----------------------------
%% Load lateral controller parameters
% ----------------------------
LLC = load_LowLevelControlData();
LLC_sampleTime = LLC.sample_time;

purePursuitParams   = purePursuitControllerParams();
stanleyParams       = stanleyControllerParams();
clothoidBasedParams = clothoidBasedControllerParams();
previewPointParams = previewPointControllerParams();

% ----------------------------
%% Select lateral controller type
% ----------------------------
select_lateralController;

%------------------------------------------------------------------------------
%% Longitudinal velocity estimation parameters
%------------------------------------------------------------------------------

omega_LPF = 100;
vmin      = 1;
delta     = 0.1;
beta      = 0.8;
hv        = 0.2;
ha        = 0.05;
hb        = 0.005;

%------------------------------------------------------------------------------
%% ABS parameters
%------------------------------------------------------------------------------

rho_off    = 0.9;
rho_on     = 1.02;
rho        = 0.97;
lambda1_th = 0.1;
lambda2_th = 0.1;
Tb_dot_max = 5e3;
v_on       = 2.5;
hvp        = 0.5;
v_stop     = 0.1;

%------------------------------------------------------------------------------
%% ESP proportional controller coefficient
%------------------------------------------------------------------------------

kp = 20;

% ----------------------------
%% Load road scenario
% ----------------------------
% Select scenario
scenario_name = 'abs_straight_road_1';  % Standard condition -> Objective stop before 78 m
% scenario_name = 'abs_straight_road_2';  % Wet asphalt -> Objective stop before 90 m
% scenario_name = 'abs_straight_road_3';   % Ice on left wheels -> Objective keep the vehicle in the road and stop before 120 m

scenario_data         = load(strcat('./Scenario/', scenario_name , '.mat'));
road_data         = [scenario_data.path.x, scenario_data.path.y, scenario_data.path.theta];
road_data_sampled = [scenario_data.path.x_sampled', scenario_data.path.y_sampled'];
road_condition = scenario_data.path.road_condition;


% Define initial conditions for the simulation
Ts = scenario_data.times.step_size;  % integration step for the simulation (fixed step)
T0 = scenario_data.times.t0;         % starting time of the simulation
Tf = scenario_data.times.tf;         % stop time of the simulation


% Desired vehicle speed
max_speed = scenario_data.vehicle_control.max_speed;

% type of longitudinal control
low_level_control = scenario_data.vehicle_control.low_level_control;

% Direct control configured
time_pedal = scenario_data.vehicle_control.time_pedal;
req_pedal = scenario_data.vehicle_control.req_pedal;

% Low_level_control
time_manoeuvre = scenario_data.vehicle_control.time_manoeuvre;
req_speed = scenario_data.vehicle_control.req_speed;

% ----------------------------
%% Define initial conditions for the simulation
% ----------------------------
X0 = loadInitialConditions(scenario_data.vehicle_control.max_speed/3.6);


% ----------------------------
%% Define graphical interface settings
% ----------------------------
% Set this flag to 1 in order to enable online plots during the simulation
enable_onlinePlots = 0;
% Set this flag to 1 in order to enable a zoomed view in online simulation plots
enable_zoom = 1;

if (enable_onlinePlots)
    % Initialize figure for online plots
    figure('Name','Road Scenario','NumberTitle','off')
    grid on
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    title('Road Scenario')
    hold on
end

%------------------------------------------------------------------------------
%% Start Simulation
%------------------------------------------------------------------------------
fprintf( 'Starting Simulation\n' )
fprintf( 'initial speed: %.0f (km/h)\n', max_speed );

tic; % start measuring time

% Run the simulation
model_sim = sim( 'Vehicle_Model_2Track_contABS_EMPTY' );

elapsed_time_simulation = toc; % stop measuring time

fprintf( 'Simulation completed\n' )
fprintf( 'The total simulation time was %.2f seconds\n', ...
         elapsed_time_simulation)

%------------------------------------------------------------------------------
%% Post-Processing
%------------------------------------------------------------------------------

dataAnalysis( model_sim, vehicle_data, Ts, scenario_name );


