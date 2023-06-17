function vehicle_data = getVehicleDataStruct()

% ----------------------------------------------------------------
%% Function purpose: define a struct containing vehicle data. 
%                    All parameters refer to the vehicle Chimera Evoluzione
% ----------------------------------------------------------------

% ----------------------------------------------------------------
%  _____               ___       _        
% |_   _|  _ _ _ ___  |   \ __ _| |_ __ _ 
%   | || || | '_/ -_) | |) / _` |  _/ _` |
%   |_| \_, |_| \___| |___/\__,_|\__\__,_|
%       |__/                              
% ----------------------------------------------------------------
tire.Vlow_long                 = 8;    % [m/s] speed threshold to use low-speed corrections in the tire longit slip and force models
tire.Vlow_lat                  = 8;    % [m/s] speed threshold to use low-speed corrections in the tire lateral slip and force models

load_MF96_tyre_data;


% ----------------------------------------------------------------
%  ___                           _            ___       _        
% / __|_  _ ____ __  ___ _ _  __(_)___ _ _   |   \ __ _| |_ __ _ 
% \__ \ || (_-< '_ \/ -_) ' \(_-< / _ \ ' \  | |) / _` |  _/ _` |
% |___/\_,_/__/ .__/\___|_||_/__/_\___/_||_| |___/\__,_|\__\__,_|
%             |_|                                                
% ----------------------------------------------------------------

% REAR SUSPENSIONS
rear_suspension.Ks_r      = ((26000)^(-1) + (100*10^3)^(-1))^(-1); % [N/m] Rear suspension+tire 
                                                                   % stiffness (2 springs in series)
rear_suspension.Cs_r      = 2125;  % [N*s/m] Rear suspension damping
rear_suspension.Cs_r_b    = 1750;  % [N*s/m] Rear suspension damping bound
rear_suspension.Cs_r_r    = 2500;  % [N*s/m] Rear suspension damping rebound
rear_suspension.Karb_r    = 0;     % [Nm/rad] anti-roll bar stiffness
rear_suspension.stroke_r  = 0.06;  % [m] maximum rear damper stroke
rear_suspension.K_es_r    = 50000; % [N/m] rear damper's end-stops stiffness
rear_suspension.C_es_r    = 2000;  % [N*s/m] rear damper's end-stops damping
rear_suspension.h_rc_r    = 0.033; % [m] rear roll center height   
rear_suspension.z__rlx_r  = 0.175; % [m] spring free length
rear_suspension.reg_fact  = 1e5;   % [1/m] regularized sign steepness factor (equal for front and rear)

% FRONT SUSPENSIONS
front_suspension.Ks_f     = ((20000)^(-1) + (100*10^3)^(-1))^(-1); % [N/m] Rear suspension+tire
                                                                   % stiffness (2 springs in series)
front_suspension.Cs_f     = 2125;         % [N*s/m] Front suspension dumping (mean for state space tuning)
front_suspension.Cs_f_b   = 1750;         % [N*s/m] Rear suspension damping bound
front_suspension.Cs_f_r   = 2500;         % [N*s/m] Rear suspension damping rebound
front_suspension.Karb_f   = 195.4*180/pi; % [Nm/rad] anti-roll bar stiffness
front_suspension.stroke_f = 0.06;         % [m] maximum front damper stroke
front_suspension.K_es_f   = 50000;        % [N/m] front damper's end-stops stiffness
front_suspension.C_es_f   = 2000;         % [N*s/m] front damper's end-stops damping
front_suspension.h_rc_f   = 0.056;        % [m] front roll center height     
front_suspension.z__rlx_f = 0.175;        % [m] spring free length

suspension.camber_gain    = 0.997836;     % [-] camber gain constant (linear fitting from suspension kinematic model)
suspension.tau_N          = 0.06;         % [s] time constant for the vertical loads dynamics

% ----------------------------------------------------------------
%   ___ _               _      ___       _        
%  / __| |_  __ _ _____(_)___ |   \ __ _| |_ __ _ 
% | (__| ' \/ _` (_-<_-< (_-< | |) / _` |  _/ _` |
%  \___|_||_\__,_/__/__/_/__/ |___/\__,_|\__\__,_|
%                                                 
% ----------------------------------------------------------------

% CHASSIS (including all the sprung mass)
% Inertia tensor for the chassis
% is =  |  is_xx   0   -is_xz |
%       |    0   is_yy    0   |
%       | -is_xz   0    is_zz |
chassis.is_xx = 0.9*85;     % [kg*m^2] chassis moment of inertia about x axis
chassis.is_yy = 0.9*85;     % [kg*m^2] chassis moment of inertia about y axis
chassis.is_zz = 1300;       % [kg*m^2] chassis moment of inertia about z axis
chassis.is_xz = 0.9*60;     % [kg*m^2] chassis product of inertia xz

% Load the new vehicle data for a formula SAE car
student_vehicle_data;

% ----------------------------------------------------------------
%  _   _                                   ___       _        
% | | | |_ _  ____ __ _ _ _  _ _ _  __ _  |   \ __ _| |_ __ _ 
% | |_| | ' \(_-< '_ \ '_| || | ' \/ _` | | |) / _` |  _/ _` |
%  \___/|_||_/__/ .__/_|  \_,_|_||_\__, | |___/\__,_|\__\__,_|
%               |_|                |___/                      
% ----------------------------------------------------------------

% UNSRPUNG BODY IS MADE OF 4 WHEELS, SUSPENSIONS, TRANSMISSION AND BRAKE MASSES

% WHEEL
% Inertia tensor of the wheel
% iwd = | iwd   0  0  |
%       |  0  iwa  0  |
%       |  0   0  iwd |

% REAR WHEELS
m_wr = new_vehicle.m_u_r / 2; % [kg] Wheel mass
w_wr = 6*25.4*10^(-3);  % [m] Wheel width
rr   = tyre_data_r.R0;  % [m] Rear Wheel Radius    
rear_wheel.Rr            = rr;                
rear_wheel.width         = w_wr;          
rear_wheel.mass          = m_wr;  
rear_wheel.iwd_r         = m_wr/12 * (3*rr^2 + w_wr^2);  % [kg*m^2] inertia of the wheel
rear_wheel.iwa_r         = 0.5;  % [kg*m^2] inertia of the whole wheel assembly
rear_wheel.static_camber = new_vehicle.gamma_r;  % [deg] Static camber for rear wheels

m_ur = 2*m_wr;                   % [kg] Rear unsprung mass 
rear_unsprung.mass = m_ur;   

% FRONT WHEELS
m_wf = new_vehicle.m_uf / 2;               % [kg] Wheel mass
w_wf = 6*25.4*10^(-3);  % [m] Wheel width
rf   = tyre_data_f.R0;  % [m] Front Wheel Radius
front_wheel.Rf            = rf;                       
front_wheel.width         = w_wf;                     
front_wheel.mass          = m_wf;                     
front_wheel.iwd_f         = m_wf/12 * (3*rf^2 + w_wf^2); % [kg*m^2] inertia of the wheel 
front_wheel.iwa_f         = 0.5; % [kg*m^2] inertia of the whole wheel assembly
front_wheel.static_camber = new_vehicle.gamma_f;   % [deg] Static camber for rear wheels

m_uf = 2*m_wf;                   % [kg] Front unsprung mass
front_unsprung.mass = m_uf;  


% ----------------------------------------------------------------
%    ___                   _ _  __   __   _    _    _       ___       _        
%   / _ \__ _____ _ _ __ _| | | \ \ / /__| |_ (_)__| |___  |   \ __ _| |_ __ _ 
%  | (_) \ V / -_) '_/ _` | | |  \ V / -_) ' \| / _| / -_) | |) / _` |  _/ _` |
%   \___/ \_/\___|_| \__,_|_|_|   \_/\___|_||_|_\__|_\___| |___/\__,_|\__\__,_|
%                                                                              
% ----------------------------------------------------------------

% VEHICLE    
vehicle.Lf    = new_vehicle.Lf;  % [m] Distance between vehicle CoM and front wheels axle      
vehicle.Lr    = new_vehicle.Lr;  % [m] Distance between vehicle CoM and front wheels axle      
vehicle.L     = vehicle.Lf + vehicle.Lr;   % [m] Vehicle wheelbase  
vehicle.hGs   = new_vehicle.h_G;   % [m] CoM vertical position  
vehicle.Wf    = new_vehicle.W_f;   % [m] Front track width     
vehicle.Wr    = new_vehicle.W_r;   % [m] Rear track width  
vehicle.lx    = 0.3;    % [m] Tire relaxation length, for longitudinal slip dynamics   
vehicle.ly    = 0.3;    % [m] Tire relaxation length, for lateral slip dynamics 
% Inertia and mass
m = new_vehicle.m_s + m_uf + m_ur; % [kg] Total mass of the vehicle + driver = 245 + 70 = 315 kg
vehicle.m     = m;    
vehicle.i_xx  = new_vehicle.Jx_s;     % [kg*m^2] Moment of inertia of the vehicle w.r.t. x axis
vehicle.i_yy  = new_vehicle.Jy_s;    % [kg*m^2] Moment of inertia of the vehicle w.r.t. y axis
vehicle.i_zz  = new_vehicle.Jz_s;   % [kg*m^2] Moment of inertia of the vehicle w.r.t. z axis
vehicle.i_xz  = new_vehicle.Jxz_s;    % [kg*m^2] Product of inertia of the vehicle 
vehicle.g     = 9.81;   % [m/s^2] acceleration due to gravity


% ----------------------------------------------------------------
%   __  __                           _    ___      ___
%  |  \/  |__ _ ______  __ _ _ _  __| |  / __|___ / __|
%  | |\/| / _` (_-<_-< / _` | ' \/ _` | | (__/ _ \ (_ |
%  |_|  |_\__,_/__/__/ \__,_|_||_\__,_|  \___\___/\___|
%
% ----------------------------------------------------------------

vehicle.ms = m - m_uf - m_ur; % [kg] Sprung Mass


% ----------------------------------------------------------------
%     _                   _                      _
%    /_\  ___ _ _ ___  __| |_  _ _ _  __ _ _ __ (_)__ ___
%   / _ \/ -_) '_/ _ \/ _` | || | ' \/ _` | '  \| / _(_-<
%  /_/ \_\___|_| \___/\__,_|\_, |_||_\__,_|_|_|_|_\__/__/
%                           |__/
% ----------------------------------------------------------------

aerodynamics.CAx  = 0.8;    % [N*s^2/m^2] Aero drag coefficient
aerodynamics.CAzf = 0.4; % [N*s^2/m^2] Aero downforce coeff at front axle
aerodynamics.CAzr = 0.5; % [N*s^2/m^2] Aero downforce coeff at rear axle


% ----------------------------------------------------------------
%   _____                       _       _            ___       _        
%  |_   _| _ __ _ _ _  ____ __ (_)_____(_)___ _ _   |   \ __ _| |_ __ _ 
%    | || '_/ _` | ' \(_-< '  \| (_-<_-< / _ \ ' \  | |) / _` |  _/ _` |
%    |_||_| \__,_|_||_/__/_|_|_|_/__/__/_\___/_||_| |___/\__,_|\__\__,_|
%                                                                       
% ----------------------------------------------------------------

transmission.tau_red = 52/15;  % [-] Transmission ratio of the gearbox
transmission.eff_red = 0.93;   % [-] Efficiency of the gearbox


% ----------------------------------------------------------------
%   ___ _               _             ___             ___       _        
%  / __| |_ ___ ___ _ _(_)_ _  __ _  / __|_  _ ___   |   \ __ _| |_ __ _ 
%  \__ \  _/ -_) -_) '_| | ' \/ _` | \__ \ || (_-<_  | |) / _` |  _/ _` |
%  |___/\__\___\___|_| |_|_||_\__, | |___/\_, /__(_) |___/\__,_|\__\__,_|
%                             |___/       |__/                           
% ----------------------------------------------------------------

steering_system.tau_D = 12; %3.67;  % [-] Steering transmission ratio (pinion-rack)
steering_system.tau_H = new_vehicle.tau_H;       % [s] Time constant for steering wheel dynamics


% Braking system                        
braking.max_brake_torque_front = 3500; %600;  % [Nm] max front braking torque that the hydraulic system can provide
braking.max_brake_torque_rear  = 3500; %600;  % [Nm] max rear braking torque that the hydraulic system can provide
braking.brakeRatio             = 0.7;  % [-] front/rear brake circuits pressure distribution
braking.totBrakeTorque         = 7000; %750;  % [Nm] max total brake torque that the braking system can develop (it is then split btween front/rear axles) 
braking.tau_br                 = 0.03; % [s] time constant for brake actuation dynamics
braking.regularSignScale       = 1;    % [rad/s] scale parameter for the regularized sign function


% ----------------------------------------------------------------
%  __  __     _               ___
% |  \/  |___| |_  ___  _ _  |   \ __ _| |_ __ _
% | |\/| / _ |  _|/ _ \| '_| | |) / _` |  _/ _` |
% |_|  |_\___/\_ |\___/|_|   |___/\__,_|\__\__,_|
%
% ----------------------------------------------------------------

% Electric motor parameters (motor model: Emrax 208, with Chimera Evoluzione powertrain)
motor.maxTorque         = 100;     % [Nm] max torque that the motor can provide
motor.speedForTorqueCut = 4800;   % [rpm] motor rotational speed at which torque is decreased a lot
motor.maxRotSpeed       = 5200;   % [rpm] max rotational speed of the motor
motor.k_torque          = 0.83;   % [-] motor torque constant
motor.I_max             = 100;    % [A] max motor current
motor.tau_mot           = 0.03;   % [s] time constant for motor actuation dynamics
motor.tau_ped           = 0.08;   % [s] time constant for pedal dynamics


% ----------------------------------------------------------------
%  ___       _   _                  ___       _
% | |)) __ _| |_| |_ ___ _ _ _  _  |   \ __ _| |_ __ _
% | |)\/ _` |  _|  _/ -_) '_| || | | |) / _` |  _/ _` |
% |___/\__,_|\__\\__\___|_|  \_, | |___/\__,_|\__\__,_|
%                            |__/ 
% ----------------------------------------------------------------

accumulator.maxPower = 175;  % [kW] max output power for the battery pack


% ----------------------------------------------------------------
% __   __   _    _    _       ___ _               _     ___       _        
% \ \ / /__| |_ (_)__| |___  / __| |_ _ _ _  _ __| |_  |   \ __ _| |_ __ _ 
%  \ V / -_) ' \| / _| / -_) \__ \  _| '_| || / _|  _| | |) / _` |  _/ _` |
%   \_/\___|_||_|_\__|_\___| |___/\__|_|  \_,_\__|\__| |___/\__,_|\__\__,_|
% 
% ----------------------------------------------------------------

% Store all sub-structures of sub-systems data in vehicle structure

vehicle_data.chassis          = chassis;

vehicle_data.aerodynamics     = aerodynamics;

vehicle_data.transmission     = transmission;

vehicle_data.steering_system  = steering_system;

vehicle_data.rear_unsprung    = rear_unsprung;
vehicle_data.front_unsprung   = front_unsprung;

vehicle_data.rear_wheel       = rear_wheel;
vehicle_data.front_wheel      = front_wheel;

vehicle_data.suspension       = suspension;
vehicle_data.rear_suspension  = rear_suspension;
vehicle_data.front_suspension = front_suspension;

vehicle_data.motor            = motor;

vehicle_data.braking          = braking;

vehicle_data.tire             = tire;

vehicle_data.tyre_data_f      = tyre_data_f;
vehicle_data.tyre_data_r      = tyre_data_r;

vehicle_data.accumulator      = accumulator;

vehicle_data.vehicle          = vehicle;

end


