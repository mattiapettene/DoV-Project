function dataAnalysis(model_sim,vehicle_data,Ts,scenario_name)

% % export images
% path = 'Plots/imageNOESP';
% figure_width_l = 800;
% figure_heigth_l = figure_width_l/1.5;
% figure_pos_x_l = 200;
% figure_pos_y_l = 200;

% ----------------------------------------------------------------
%% Post-Processing and Data Analysis
% ----------------------------------------------------------------
path = load(strcat('./Scenario/', scenario_name, '.mat'));
xRoadLeft = path.path.x_left_sampled;
xRoadRight = path.path.x_right_sampled;
yRoadLeft = path.path.y_left_sampled;
yRoadRight = path.path.y_right_sampled;
xRoad = path.path.x_sampled;
yRoad = path.path.y_sampled;

% ---------------------------------
%% Load vehicle data
% ---------------------------------
Lf = vehicle_data.vehicle.Lf;  % [m] Distance between vehicle CoG and front wheels axle
Lr = vehicle_data.vehicle.Lr;  % [m] Distance between vehicle CoG and front wheels axle
L  = vehicle_data.vehicle.L;   % [m] Vehicle length
Wf = vehicle_data.vehicle.Wf;  % [m] Width of front wheels axle
Wr = vehicle_data.vehicle.Wr;  % [m] Width of rear wheels axle
m  = vehicle_data.vehicle.m;   % [kg] Vehicle Mass
g  = vehicle_data.vehicle.g;   % [m/s^2] Gravitational acceleration
tau_D = vehicle_data.steering_system.tau_D;  % [-] steering system ratio (pinion-rack)
Rf_wheel = vehicle_data.front_wheel.Rf; %[m] Front wheel radius
Rr_wheel = vehicle_data.rear_wheel.Rr; %[m] Rear wheel radius


% ---------------------------------
%% Extract data from simulink model
% ---------------------------------
time_sim = model_sim.states.u.time;
dt = time_sim(2)-time_sim(1);

% -----------------
% Inputs
% -----------------
ped_0      = model_sim.inputs.ped_0.data;
delta_D    = model_sim.inputs.delta_D.data;

% -----------------
% States
% -----------------
x_CoM      = model_sim.states.x.data;
y_CoM      = model_sim.states.y.data;
psi        = model_sim.states.psi.data;
u          = model_sim.states.u.data;
v          = model_sim.states.v.data;
Omega      = model_sim.states.Omega.data;
Fz_rr      = model_sim.states.Fz_rr.data;
Fz_rl      = model_sim.states.Fz_rl.data;
Fz_fr      = model_sim.states.Fz_fr.data;
Fz_fl      = model_sim.states.Fz_fl.data;
delta      = model_sim.states.delta.data;
omega_rr   = model_sim.states.omega_rr.data;
omega_rl   = model_sim.states.omega_rl.data;
omega_fr   = model_sim.states.omega_fr.data;
omega_fl   = model_sim.states.omega_fl.data;
alpha_rr   = model_sim.states.alpha_rr.data;
alpha_rl   = model_sim.states.alpha_rl.data;
alpha_fr   = model_sim.states.alpha_fr.data;
alpha_fl   = model_sim.states.alpha_fl.data;
kappa_rr   = model_sim.states.kappa_rr.data;
kappa_rl   = model_sim.states.kappa_rl.data;
kappa_fr   = model_sim.states.kappa_fr.data;
kappa_fl   = model_sim.states.kappa_fl.data;

% -----------------
% Extra Parameters
% -----------------
Tw_rr      = model_sim.extra_params.Tw_rr.data;
Tw_rl      = model_sim.extra_params.Tw_rl.data;
Tw_fr      = model_sim.extra_params.Tw_fr.data;
Tw_fl      = model_sim.extra_params.Tw_fl.data;
Fx_rr      = model_sim.extra_params.Fx_rr.data;
Fx_rl      = model_sim.extra_params.Fx_rl.data;
Fx_fr      = model_sim.extra_params.Fx_fr.data;
Fx_fl      = model_sim.extra_params.Fx_fl.data;
Fy_rr      = model_sim.extra_params.Fy_rr.data;
Fy_rl      = model_sim.extra_params.Fy_rl.data;
Fy_fr      = model_sim.extra_params.Fy_fr.data;
Fy_fl      = model_sim.extra_params.Fy_fl.data;
Mz_rr      = model_sim.extra_params.Mz_rr.data;
Mz_rl      = model_sim.extra_params.Mz_rl.data;
Mz_fr      = model_sim.extra_params.Mz_fr.data;
Mz_fl      = model_sim.extra_params.Mz_fl.data;
gamma_rr   = model_sim.extra_params.gamma_rr.data;
gamma_rl   = model_sim.extra_params.gamma_rl.data;
gamma_fr   = model_sim.extra_params.gamma_fr.data;
gamma_fl   = model_sim.extra_params.gamma_fl.data;
delta_fr   = model_sim.extra_params.delta_fr.data;
delta_fl   = model_sim.extra_params.delta_fl.data;



% Chassis side slip angle beta [rad]
beta = atan(v./u);

% -----------------
% Accelerations
% -----------------
% Derivatives of u, v [m/s^2]
dot_u = diff(u)/Ts;
dot_v = diff(v)/Ts;
% Total longitudinal and lateral accelerations
Ax = dot_u(1:end) - Omega(2:end).*v(2:end);
Ay = dot_v(1:end) + Omega(2:end).*u(2:end);
% Ax low-pass filtered signal (zero-phase digital low-pass filtering)
Wn_filter = 0.01;
[b_butt,a_butt] = butter(4,Wn_filter,'low');
Ax_filt = filtfilt(b_butt,a_butt,Ax);
dot_u_filt = filtfilt(b_butt,a_butt,dot_u);
% Steady state lateral acceleration
Ay_ss = Omega.*u;
% Longitudinal jerk [m/s^3]
jerk_x = diff(dot_u)/Ts;

% -----------------
% Other parameters
% -----------------
% Total CoM speed [m/s]
vG = sqrt(u.^2 + v.^2);
% Steady state and transient curvature [m]
rho_ss   = Omega./vG;
rho_tran = ((dot_v.*u(1:end-1) - dot_u.*v(1:end-1)) ./ ((vG(1:end-1)).^3)) + rho_ss(1:end-1);
% Desired sinusoidal steering angle for the equivalent single track front wheel
desired_steer_atWheel = delta_D/tau_D;


% ---------------------------------
%% PLOTS
% ---------------------------------

% ---------------------------------
%% Plot vehicle inputs
% ---------------------------------
figure('Name','Inputs','NumberTitle','off'), clf
% --- pedal --- %
ax(1) = subplot(211);
hold on
plot(time_sim,ped_0,'LineWidth',2)
grid on
title('pedal $p_0$ [-]')
xlim([0 time_sim(end)])
% --- delta_0 --- %
ax(2) = subplot(212);
plot(time_sim,delta_D,'LineWidth',2)
grid on
title('steering angle $\delta_D$ [deg]')
xlim([0 time_sim(end)])

% ---------------------------------
%% Plot vehicle motion
% ---------------------------------
figure('Name','Vehicle otion','NumberTitle','off'), clf
% --- u --- %
ax(1) = subplot(221);
plot(time_sim,u*3.6,'LineWidth',2)
grid on
title('$u$ [km/h]')
xlim([0 time_sim(end)])
% --- v --- %
ax(2) = subplot(222);
plot(time_sim,v,'LineWidth',2)
grid on
title('$v$ [m/s]')
xlim([0 time_sim(end)])
% --- Omega --- %
ax(3) = subplot(223);
plot(time_sim,Omega,'LineWidth',2)
grid on
title('$\Omega$ [rad/s]')
xlim([0 time_sim(end)])

% ---------------------------------
%% Plot steering angles
% ---------------------------------
figure('Name','Steer','NumberTitle','off'), clf
% --- delta_0 --- %
ax(1) = subplot(211);
plot(time_sim,delta_D,'LineWidth',2)
grid on
title('$\delta_0$ [deg]')
xlim([0 time_sim(end)])
% --- comparison --- %
ax(4) = subplot(212);
hold on
title('$\delta$ [deg]')
plot(time_sim,delta_D/tau_D,'LineWidth',2)
plot(time_sim,delta_fr,'LineWidth',2)
plot(time_sim,delta_fl,'LineWidth',2)
grid on
legend('$\delta_D/\tau_D$','$\delta_{fr}$','$\delta_{fl}$','location','best')
xlim([0 time_sim(end)])

% -------------------------------
%% Plot lateral tire slips and lateral forces
% -------------------------------
figure('Name','Lateral slips & forces','NumberTitle','off'), clf
% --- alpha --- %
ax(1) = subplot(211);
hold on
plot(time_sim,alpha_rr,'LineWidth',2)
plot(time_sim,alpha_rl,'LineWidth',2)
plot(time_sim,alpha_fr,'LineWidth',2)
plot(time_sim,alpha_fl,'LineWidth',2)
grid on
title('$\alpha$ [deg]')
legend('$\alpha_{rr}$ [deg]','$\alpha_{rl}$ [deg]','$\alpha_{fr}$ [deg]','$\alpha_{fl}$ [deg]');
xlim([0 time_sim(end)])
% --- Fy --- %
ax(2) = subplot(212);
hold on;
plot(time_sim,Fy_rr,'LineWidth',2)
plot(time_sim,Fy_rl,'LineWidth',2)
plot(time_sim,Fy_fr,'LineWidth',2)
plot(time_sim,Fy_fl,'LineWidth',2)
grid on
title('$Fy$ [N]')
legend('$Fy_{rr}$ [N]','$Fy_{rl}$ [Nm]','$Fy_{fr}$ [N]','$Fy_{fl}$ [N]')
xlim([0 time_sim(end)])
% linkaxes(ax,'x')
clear ax

% ---------------------------------
%% Plot longitudinal tire slips and longitudinal forces
% ---------------------------------
figure('Name','Long slips & forces','NumberTitle','off'), clf
% --- kappa --- %
ax(1) = subplot(211);
hold on;
plot(time_sim,kappa_rr,'LineWidth',2)
plot(time_sim,kappa_rl,'LineWidth',2)
plot(time_sim,kappa_fr,'LineWidth',2)
plot(time_sim,kappa_fl,'LineWidth',2)
grid on
title('$\kappa$ [-]')
legend('$\kappa_{rr}$ [-]','$\kappa_{rl}$ [-]','$\kappa_{fr}$ [-]','$\kappa_{fl}$ [-]')
xlim([0 time_sim(end)])

% --- Fx --- %
ax(2) = subplot(212);
hold on;
plot(time_sim,Fx_rr,'LineWidth',2)
plot(time_sim,Fx_rl,'LineWidth',2)
plot(time_sim,Fx_fr,'LineWidth',2)
plot(time_sim,Fx_fl,'LineWidth',2)
grid on
title('$Fx$ [N]')
legend('$Fx_{rr}$ [N]','$Fx_{rl}$ [N]','$Fx_{fr}$ [N]','$Fx_{fl}$ [N]')
xlim([0 time_sim(end)])
clear ax

% hgexport(f, sprintf('%s/LongitudinalSlip.eps', path))

% ---------------------------------
%% Plot wheel torques and wheel rates
% ---------------------------------
figure('Name','Wheel rates & torques','NumberTitle','off'), clf
% --- omega --- %
ax(1) = subplot(211);
hold on
plot(time_sim,omega_rr,'LineWidth',2)
plot(time_sim,omega_rl,'LineWidth',2)
plot(time_sim,omega_fr,'LineWidth',2)
plot(time_sim,omega_fl,'LineWidth',2)
grid on
title('$\omega$ [rad/s]')
legend('$\omega_{rr}$ [rad/s]','$\omega_{rl}$ [rad/s]','$\omega_{fr}$ [rad/s]','$\omega_{fl}$ [rad/s]')
xlim([0 time_sim(end)])

% --- Tw --- %
ax(2) = subplot(212);
hold on
plot(time_sim,Tw_rr,'LineWidth',2)
plot(time_sim,Tw_rl,'LineWidth',2)
plot(time_sim,Tw_fr,'LineWidth',2)
plot(time_sim,Tw_fl,'LineWidth',2)
grid on
title('$Tw$ [Nm]')
legend('$Tw_{rr}$ [Nm]','$Tw_{rl}$ [Nm]','$Tw_{fr}$ [Nm]','$Tw_{fl}$ [Nm]')
xlim([0 time_sim(end)])

% linkaxes(ax,'x')
clear ax

% ---------------------------------
%% Plot vertical tire loads and self-aligning torques
% ---------------------------------
figure('Name','Vert loads & aligning torques','NumberTitle','off'), clf
% --- Fz --- %
ax(1) = subplot(211);
hold on
plot(time_sim,Fz_rr,'LineWidth',2)
plot(time_sim,Fz_rl,'LineWidth',2)
plot(time_sim,Fz_fr,'LineWidth',2)
plot(time_sim,Fz_fl,'LineWidth',2)
grid on
title('$Fz_{rr}$ [N]')
legend('$Fz_{rr}$ [N]','$Fz_{rl}$ [N]','$Fz_{fr}$ [N]','$Fz_{fl}$ [N]')
xlim([0 time_sim(end)])
% --- Mz --- %
ax(2) = subplot(212);
hold on
plot(time_sim,Mz_rr,'LineWidth',2)
plot(time_sim,Mz_rl,'LineWidth',2)
plot(time_sim,Mz_fr,'LineWidth',2)
plot(time_sim,Mz_fl,'LineWidth',2)
grid on
title('$Mz$ [Nm]')
legend('$Mz_{rr}$ [Nm]','$Mz_{rl}$ [Nm]','$Mz_{fr}$ [Nm]','$Mz_{fl}$ [Nm]')
xlim([0 time_sim(end)])

% linkaxes(ax,'x')
clear ax

% ---------------------------------
%% Plot wheel camber
% ---------------------------------
figure('Name','Wheel camber','NumberTitle','off'), clf
% --- gamma_rr --- %
ax(1) = subplot(221);
plot(time_sim,gamma_rr,'LineWidth',2)
grid on
title('$\gamma_{rr}$ [deg]')
xlim([0 time_sim(end)])
% --- gamma_rl --- %
ax(2) = subplot(222);
plot(time_sim,gamma_rl,'LineWidth',2)
grid on
title('$\gamma_{rl}$ [deg]')
xlim([0 time_sim(end)])
% --- gamma_fr --- %
ax(3) = subplot(223);
plot(time_sim,gamma_fr,'LineWidth',2)
grid on
title('$\gamma_{fr}$ [deg]')
xlim([0 time_sim(end)])
% --- gamma_fl --- %
ax(4) = subplot(224);
plot(time_sim,gamma_fl,'LineWidth',2)
grid on
title('$\gamma_{fl}$ [deg]')
xlim([0 time_sim(end)])

% linkaxes(ax,'x')
clear ax

% ---------------------------------
%% Plot accelerations, chassis side slip angle and curvature
% ---------------------------------
figure('Name','Pars extra','NumberTitle','off'), clf
% --- ax --- %
ax(1) = subplot(221);
plot(time_sim(2:end),dot_u - Omega(2:end).*v(2:end),'LineWidth',2)
hold on
plot(time_sim(2:end),diff(u)/Ts,'--g','LineWidth',2)
plot(time_sim(2:end),Ax_filt,'-.b','LineWidth',1)
plot(time_sim(2:end),dot_u_filt,'-.r','LineWidth',1)
grid on
title('$a_{x}$ $[m/s^2]$')
legend('$\dot{u}-\Omega v$','$\dot{u}$','filt $\dot{u}-\Omega v$','filt $\dot{u}$','Location','northeast')
xlim([0 time_sim(end)])
% --- ay --- %
ax(2) = subplot(222);
plot(time_sim(2:end),dot_v + Omega(2:end).*u(2:end),'LineWidth',2)
hold on
plot(time_sim(2:end),Omega(2:end).*u(2:end),'--g','LineWidth',1)
grid on
title('$a_{y}$ $[m/s^2]$')
legend('$\dot{v}+\Omega u$','$\Omega u$','Location','best')
xlim([0 time_sim(end)])
% --- beta --- %
ax(3) = subplot(223);
plot(time_sim,rad2deg(beta),'LineWidth',2)
grid on
title('$\beta$ [deg]')
xlim([0 time_sim(end)])
% --- rho --- %
ax(4) = subplot(224);
plot(time_sim,rho_ss,'LineWidth',2)
hold on
plot(time_sim(1:end-1),rho_tran,'--g','LineWidth',1)
grid on
title('$\rho$ [$m^{-1}$]')
legend('$\rho_{ss}$','$\rho_{transient}$','Location','best')
xlim([0 time_sim(end)])

% linkaxes(ax,'x')
clear ax

% ---------------------------------
%% Plot vehicle pose x,y,psi
% ---------------------------------
f = figure('Name','Pose','NumberTitle','off');
% f.Position = [figure_pos_x_l figure_pos_y_l figure_width_l figure_heigth_l];
% --- x --- %
ax(1) = subplot(221);
plot(time_sim,x_CoM,'LineWidth',2)
grid on
title('$x$ [m]')
xlim([0 time_sim(end)])
% --- y --- %
ax(2) = subplot(222);
plot(time_sim,y_CoM,'LineWidth',2)
grid on
title('$y$ [m]')
xlim([0 time_sim(end)])
% --- psi --- %
ax(3) = subplot(223);
plot(time_sim,rad2deg(psi),'LineWidth',2)
grid on
title('$\psi$ [deg]')
xlim([0 time_sim(end)])
% hgexport(f, sprintf('%s/Pose.eps', path))

% linkaxes(ax,'x')
clear ax

% -------------------------------
%% Plot G-G diagram from simulation data
% -------------------------------
f = figure('Name','G-G plot','NumberTitle','off');
% f.Position = [figure_pos_x_l figure_pos_y_l figure_width_l figure_heigth_l];
axis equal
hold on
plot3(Ay,Ax_filt,u(1:end-1),'Color',color('purple'),'LineWidth',3)
xlabel('$a_y$ [m/s$^2$]')
ylabel('$a_x$ [m/s$^2$]')
zlabel('$u$ [m/s]')
title('G-G diagram from simulation data','FontSize',18)
grid on
% hgexport(f, sprintf('%s/GGDiagram.eps', path))

% -------------------------------
%% Plot vehicle velocity
% -------------------------------
figure('Name','Wheel rates & vehicle vel','NumberTitle','off'), clf
% --- vel --- %
hold on
plot(time_sim,u,'LineWidth',3)
plot(time_sim,omega_rr*Rr_wheel,'LineWidth',2)
plot(time_sim,omega_rl*Rr_wheel,'LineWidth',2)
plot(time_sim,omega_fr*Rf_wheel,'LineWidth',2)
plot(time_sim,omega_fl*Rf_wheel,'LineWidth',2)
grid on
title('Vehicle velocity [m/s]')
legend('$u$ [m/s]','$\omega_{rr}*Rr$ [m/s]','$\omega_{rl}*Rr$ [m/s]','$\omega_{fr}*Rf$ [m/s]','$\omega_{fl}*Rf$ [m/s]')
xlim([0 time_sim(end)])

% -------------------------------
%% Plot vehicle path
% -------------------------------
N = length(time_sim);
f = figure('Name','Real Vehicle Path','NumberTitle','off');
% f.Position = [figure_pos_x_l figure_pos_y_l figure_width_l figure_heigth_l];
set(gca,'fontsize',16)
hold on
axis equal
xlabel('x-coord [m]')
ylabel('y-coord [m]')
title('Real Vehicle Path','FontSize',18)
plot(x_CoM,y_CoM,'Color',color('gold'),'LineWidth',2)

plot(xRoadLeft,yRoadLeft,'Color','black','LineWidth',2);
plot(xRoad,yRoad,'Color','black','LineWidth',2, 'LineStyle', '--');
plot(xRoadRight,yRoadRight,'Color','black','LineWidth',2);

for i = 1:floor(N/20):N
    rot_mat = [cos(psi(i)) -sin(psi(i)) ; sin(psi(i)) cos(psi(i))];
    pos_rr = rot_mat*[-Lr -Wr/2]';
    pos_rl = rot_mat*[-Lr +Wr/2]';
    pos_fr = rot_mat*[+Lf -Wf/2]';
    pos_fl = rot_mat*[+Lf +Wf/2]';
    pos = [pos_rr pos_rl pos_fl pos_fr];
    p = patch(x_CoM(i) + pos(1,:),y_CoM(i) + pos(2,:),'blue');
    quiver(x_CoM(i), y_CoM(i), u(i)*cos(psi(i)), u(i)*sin(psi(i)), 'color', [1,0,0]);
    quiver(x_CoM(i), y_CoM(i), -v(i)*sin(psi(i)), v(i)*cos(psi(i)), 'color', [0.23,0.37,0.17]);
end
grid on
hold off
% hgexport(f, sprintf('%s/Path.eps', path))


% -------------------------------
%% Print Braking Distance
% -------------------------------
%brake_start = x_CoM(time_sim == 1);
%brake_stop = x_CoM(end);
fprintf('Travel Distance: %.2f (m)\n', x_CoM(end));

end

