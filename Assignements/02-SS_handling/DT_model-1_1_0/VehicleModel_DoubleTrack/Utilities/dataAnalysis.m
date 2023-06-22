function dataAnalysis(model_sim,vehicle_data,Ts)

    % ----------------------------------------------------------------
    %% Post-Processing and Data Analysis
    % ----------------------------------------------------------------

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


    %---------------------------------
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
    figure('Name','veh motion','NumberTitle','off'), clf   
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
    figure('Name','steer','NumberTitle','off'), clf   
    % --- delta_0 --- %
    ax(1) = subplot(221);
    plot(time_sim,delta_D,'LineWidth',2)
    grid on
    title('$\delta_0$ [deg]')
    xlim([0 time_sim(end)])
    % --- delta_fr --- %
    ax(2) = subplot(222);
    plot(time_sim,delta_fr,'LineWidth',2)
    grid on
    title('$\delta_{fr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- delta_fl --- %
    ax(3) = subplot(223);
    hold on
    plot(time_sim,delta_fl,'LineWidth',2)
    grid on
    title('$\delta_{fl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- comparison --- %
    ax(4) = subplot(224);
    hold on
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
    % --- alpha_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,alpha_rr,'LineWidth',2)
    grid on
    title('$\alpha_{rr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,alpha_rl,'LineWidth',2)
    grid on
    title('$\alpha_{rl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,alpha_fr,'LineWidth',2)
    grid on
    title('$\alpha_{fr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,alpha_fl,'LineWidth',2)
    grid on
    title('$\alpha_{fl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- Fy_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Fy_rr,'LineWidth',2)
    grid on
    title('$Fy_{rr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fy_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Fy_rl,'LineWidth',2)
    grid on
    title('$Fy_{rl}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Fy_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Fy_fr,'LineWidth',2)
    grid on
    title('$Fy_{fr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fy_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Fy_fl,'LineWidth',2)
    grid on
    title('$Fy_{fl}$ [N]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    % ---------------------------------
    %% Plot longitudinal tire slips and longitudinal forces
    % ---------------------------------
    figure('Name','Long slips & forces','NumberTitle','off'), clf
    % --- kappa_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,kappa_rr,'LineWidth',2)
    grid on
    title('$\kappa_{rr}$ [-]')
    xlim([0 time_sim(end)])
    % --- kappa_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,kappa_rl,'LineWidth',2)
    grid on
    title('$\kappa_{rl}$ [-]')
    xlim([0 time_sim(end)])
    % --- kappa_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,kappa_fr,'LineWidth',2)
    grid on
    title('$\kappa_{fr}$ [-]')
    xlim([0 time_sim(end)])
    % --- kappa_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,kappa_fl,'LineWidth',2)
    grid on
    title('$\kappa_{fl}$ [-]')
    xlim([0 time_sim(end)])
    % --- Fx_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Fx_rr,'LineWidth',2)
    grid on
    title('$Fx_{rr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fx_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Fx_rl,'LineWidth',2)
    grid on
    title('$Fx_{rl}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fx_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Fx_fr,'LineWidth',2)
    grid on
    title('$Fx_{fr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fx_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Fx_fl,'LineWidth',2)
    grid on
    title('$Fx_{fl}$ [N]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    % ---------------------------------
    %% Plot wheel torques and wheel rates
    % ---------------------------------
    figure('Name','Wheel rates & torques','NumberTitle','off'), clf
    % --- omega_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,omega_rr,'LineWidth',2)
    grid on
    title('$\omega_{rr}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- omega_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,omega_rl,'LineWidth',2)
    grid on
    title('$\omega_{rl}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- omega_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,omega_fr,'LineWidth',2)
    grid on
    title('$\omega_{fr}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- omega_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,omega_fl,'LineWidth',2)
    grid on
    title('$\omega_{fl}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- Tw_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Tw_rr,'LineWidth',2)
    grid on
    title('$Tw_{rr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Tw_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Tw_rl,'LineWidth',2)
    grid on
    title('$Tw_{rl}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Tw_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Tw_fr,'LineWidth',2)
    grid on
    title('$Tw_{fr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Tw_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Tw_fl,'LineWidth',2)
    grid on
    title('$Tw_{fl}$ [Nm]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    % ---------------------------------
    %% Plot vertical tire loads and self-aligning torques
    % ---------------------------------
    figure('Name','Vert loads & aligning torques','NumberTitle','off'), clf
    % --- Fz_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,Fz_rr,'LineWidth',2)
    grid on
    title('$Fz_{rr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fz_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,Fz_rl,'LineWidth',2)
    grid on
    title('$Fz_{rl}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fz_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,Fz_fr,'LineWidth',2)
    grid on
    title('$Fz_{fr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fz_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,Fz_fl,'LineWidth',2)
    grid on
    title('$Fz_{fl}$ [N]')
    xlim([0 time_sim(end)])
    % --- Mz_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Mz_rr,'LineWidth',2)
    grid on
    title('$Mz_{rr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Mz_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Mz_rl,'LineWidth',2)
    grid on
    title('$Mz_{rl}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Mz_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Mz_fr,'LineWidth',2)
    grid on
    title('$Mz_{fr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Mz_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Mz_fl,'LineWidth',2)
    grid on
    title('$Mz_{fl}$ [Nm]')
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
    figure('Name','Pose','NumberTitle','off'), clf 
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

    % linkaxes(ax,'x')
    clear ax

    % -------------------------------
    %% Plot G-G diagram from simulation data
    % -------------------------------
    figure('Name','G-G plot','NumberTitle','off'), clf
    axis equal
    hold on
    plot3(Ay,Ax_filt,u(1:end-1),'Color',color('purple'),'LineWidth',3)
    xlabel('$a_y$ [m/s$^2$]')
    ylabel('$a_x$ [m/s$^2$]')
    zlabel('$u$ [m/s]')
    title('G-G diagram from simulation data','FontSize',18)
    grid on

    % -------------------------------
    %% Plot vehicle path
    % -------------------------------
    N = length(time_sim);
    figure('Name','Real Vehicle Path','NumberTitle','off'), clf
    set(gca,'fontsize',16)
    hold on
    axis equal
    xlabel('x-coord [m]')
    ylabel('y-coord [m]')
    title('Real Vehicle Path','FontSize',18)
    plot(x_CoM,y_CoM,'Color',color('gold'),'LineWidth',2)
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
    
      %% Plot lateral load transfer
    % -------------------------------
    % Expressions used to obtain the lateral load trasfer at front and rear

    mr = m*Lf/L;
    mf = m*Lr/L;

    hrr =  vehicle_data.rear_suspension.h_rc_r; % rolling arm rear
    hrf = vehicle_data.front_suspension.h_rc_f; % rolling arm front
    hGs = vehicle_data.vehicle.hGs; % distance roll axis-com

    hr = hrf + (hrr-hrf)*Lf/(Lr+Lf);  % [m] height from ground of the projection of the CoM G on the roll axis
    hs = hGs - hr; %Distance of the rolling axis from the ground (of the com)


    % Normalized stiffness of suspension
    epsilon_phi = (vehicle_data.front_suspension.Ks_f )/((vehicle_data.front_suspension.Ks_f ) + (vehicle_data.rear_suspension.Ks_r));

    
    % RIGID CHASSIS DOUBLE TRACK-SOLUTION

    delta_f_z_r = m * Ay * ((hs*(1-epsilon_phi)/Wr) + (Lf * hrr)/(Wr * L));
    delta_f_z_r_ss = m * Ay_ss * ((hs*(1-epsilon_phi)/Wr) + (Lf * hrr)/(Wr * L));
    
    delta_f_z_f = m * Ay * ((hs *epsilon_phi)/(Wf) + (Lr * hrf)/(Wf * L));
    delta_f_z_f_ss = m * Ay_ss * ((hs *epsilon_phi)/(Wf) + (Lr * hrf)/(Wf * L));


%   Plot the nominal  -------------------------------
    figure('Name','Lateral load transfer plot','NumberTitle','off'), clf
    hold on
    plot(time_sim(2:end), delta_f_z_r,'Color',color('red'),'LineWidth',3)
    xlabel('$t$ [s]')
    ylabel('$delta_fz$ [N]')
    title('Lateral load transfer plot','FontSize',18)
    grid on
    plot(time_sim(2:end), delta_f_z_f,'Color',color('blue'),'LineWidth',3)

    plot(time_sim(2:end), (Fz_rr(2:end) - Fz_rl(2:end))/2 ,'--', 'Color','black', 'LineWidth',2)
    plot(time_sim(2:end), (Fz_fr(2:end) - Fz_fl(2:end))/2 ,'--', 'Color', 'green', 'LineWidth',2)
    hold off
    legend('Deltafzr theoretical', 'Deltafzf theoretical', 'Deltafzr measured', 'Deltafzf measured','Location', 'southeast')

    % Plot the nominal IN FUNCTION OF THE LATERAL ACCELERATION  -------------------------------
    figure('Name','Lateral load transfer plot ACCELERATION ON X','NumberTitle','off'), clf

    hold on
    plot(Ay, delta_f_z_r,'Color',color('red'),'LineWidth',3)
    xlabel('$ay [m/s^2]$')
    ylabel('$\delta_fz$ [N]')
    title('Lateral load transfer plot','FontSize',18)
    grid on
    plot(Ay, delta_f_z_f,'Color',color('blue'),'LineWidth',3)

    plot(Ay, (Fz_rr(2:end) - Fz_rl(2:end))/2 ,'--', 'Color','black', 'LineWidth',2)
    plot(Ay, (Fz_fr(2:end) - Fz_fl(2:end))/2 ,'--', 'Color', 'green', 'LineWidth',2)
    hold off
    legend('Deltafzr theoretical', 'Deltafzf theoretical', 'Deltafzr measured', 'Deltafzf measured','Location', 'southeast')

% Plot vertical tire loads difference 
    % ---------------------------------
    
    figure('Name','Veridications delta fz','NumberTitle','off'), clf
    % --- Fz_rr --- %
    subplot(1, 2, 1);
    plot(time_sim(2:end), Fz_rr(2:end) - Fz_rl(2:end) -2*delta_f_z_r ,'LineWidth',2)
    grid on
    title('$Verification Fz_{r}$ [N]')
    xlabel('$t$ [s]')
    xlim([0 time_sim(end)])
    % --- Fz_rl --- %
    
    subplot(1, 2, 2);
    plot(time_sim(2:end), Fz_fr(2:end) - Fz_fl(2:end) -2*delta_f_z_f ,'LineWidth',2)
    grid on
    title('$Verification Fz_{f}$ [N]')
    xlabel('$t$ [s]')
    xlim([0 time_sim(end)])
    
    % ---------------------------------
    
%% Normalized axle characteristic


    % Side slips - from double track
    alphaR_dt = deg2rad((alpha_rr + alpha_rl)/2);
    alphaF_dt = deg2rad((alpha_fr + alpha_fl)/2);
    delta_alpha_dt = alphaR_dt - alphaF_dt;

    % Side slips - single track
    delta_st = deg2rad((delta_fr + delta_fl)/2);
    alphaR_st = - (beta) + (Omega./u) * Lr;
    alphaF_st = delta_st - (beta) - (Omega./u) * Lf;
    delta_alpha_st = alphaR_st - alphaF_st;


    figure('Name','Side slips single track','NumberTitle','off'), clf
    hold on

    plot(time_sim, alphaR_dt, 'LineWidth',2)
    xlim([0 time_sim(end)])

    plot(time_sim, alphaF_dt, 'LineWidth',2)
    xlim([0 time_sim(end)])

    plot(time_sim, alphaR_st, '--', 'LineWidth', 1)
    xlim([0 time_sim(end)])

    plot(time_sim, alphaF_st, '--', 'LineWidth', 1)
    xlim([0 time_sim(end)])

    grid on
    legend({'$\alpha_{R}$ double track','$\alpha_{F}$ double track', '$\alpha_{R}$ single track', '$\alpha_{F}$ single track'})
    xlabel('$t$ [s]')
    ylabel('$\alpha_{R}$, $\alpha_{F}$ [rad]')

    title('Side slips $\alpha_{R}, \alpha_{F}$')

    % ------------------------------------------------------------------
    
    % Lateral forces - from double track
    
    % Real vertical load 
    Fz0R = Fz_rr + Fz_rl;
    Fz0F = Fz_fr + Fz_fl;
    
    % Static vetical load
    Fzr_0_static = m * (Lf/L) * g;
    Fzf_0_static = m * (Lr/L) * g;

    % Real forces on the axis
    Fyr_dt = Fy_rl + Fy_rr;
    Fyf_dt = sin(delta_fl).*Fx_fl + Fy_fl + sin(delta_fr).*Fx_fr + Fy_fr;
    
    % Static load normalization
    Fyr_dt_norm = Fyr_dt./Fzr_0_static;
    Fyf_dt_norm = Fyf_dt./Fzf_0_static;


    figure('Name','Normalized lateral forces','NumberTitle','off'), clf
    hold on

    plot(alphaR_dt,Fyr_dt_norm,'LineWidth',2)

    plot(alphaF_dt,Fyf_dt_norm,'LineWidth',2)

    grid on
    legend({'$Fyr$','$Fyf$'})
    xlabel('$\alpha_{R}$, $\alpha_{F}$ [rad]')
    ylabel('$Fyr/Fz0$, $Fyf/Fz0$ [-]')

    title('Normalized lateral forces')

    % ------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------
    % Handling diagram 

    % Compute the difference DeltaAlpha between rear and front side slip
    % angle
    % figure('Name','Handling diagram [Rad] ','NumberTitle','off'), clf
    % 
    % % size(Ay)
    % % size(delta_alpha_dt)
    % 
    % x1 = Ay(80000)/g;
    % y1 = -delta_alpha_dt(79999);
    % 
    % x2 = Ay(80001)/g;
    % y2 = -delta_alpha_dt(80000);
    % 
    % coefficients = polyfit([x1, x2], [y1, y2], 1);
    % slope = coefficients(1);
    % fprintf('the slope is %f rad NEW \n', slope);
    % intercept = coefficients(2);
    % 
    % % Creazione della retta
    % y = slope * (Ay/g) + intercept;
    % 
    % % Plot dei punti e della retta
    % plot(Ay/g, y, 'r');  % Retta
    % hold on;
    % scatter([x1, x2], [y1, y2], 'red', 'filled');  % Punti
    % 
    % 
    % plot(Ay/g, -delta_alpha_dt(2:end),'Color',color('blue'),'LineWidth',2)
    % 
    % plot(Ay/g, y, 'Color',color('green'),'LineWidth',2);
    % title('Handling diagram')
    % ylabel('$-Delta Alpha$ [rad]')
    % xlabel('$Ay/g$ [-]')
    % grid on
    % legend({'$-DeltaAlpha$','$tangent$'})
    % hold off;




    %-------------------------------------------------------------------
    %% Handling diagram 
    figure('Name','Handling diagram [rad] FITTING ','NumberTitle','off'), clf

    % Cut vectors
    cut_value_start = 0.05; %Selection of the starting linearizing point (Normalized acceleration value)
    index_start = find((Ay/g) > cut_value_start);
    cut_index_start = index_start(1) - 1; %Selection of the starting linearizing point (index value)

    cut_value_end = 0.4; %Selection of the ending linearizing point (Normalized acceleration value) --> THIS THE a_linear_lim
    index_end = find((Ay/g) > cut_value_end);
    cut_index_end = index_end(1) - 1; %Selection of the ending linearizing point (index value)

  
    fprintf('Il cut index  start vale %d \n', cut_index_start);
    fprintf('Il cut index  end vale %d \n', cut_index_end);
   

    % Fitting to compute the tangent (LINEAR ZONE)
    x_cut_l = Ay(cut_index_start:cut_index_end)/g;
    y_cut_l = -delta_alpha_dt((cut_index_start-1):(cut_index_end-1));

    coefficients = polyfit(x_cut_l, y_cut_l, 1);
    
    fprintf('Il coefficiente 1 della parte lineare vale = %f\n', coefficients(1));
    fprintf('Il coefficiente 0 della parte lineare vale = %f\n', coefficients(2));
    
        
    slope = coefficients(1);
    intercept = coefficients(2);

    % Creazione della retta
    y = slope * (Ay/g) + intercept;
    fprintf('Il kus calcolato nella regione lineare del fitting vale %f\n', slope);
    plot(Ay/g, zeros(size(Ay)),'Color', color('red'),'LineWidth',2);
    plot(Ay/g, -delta_alpha_dt(2:end),'Color',color('blue'),'LineWidth',2);
    hold on;

    plot(Ay/g, y, 'Color',color('green'),'LineWidth',2);
    title('Handling diagram')
    ylabel('$-Delta Alpha$ [rad]')
    xlabel('$Ay/g$ [-]')
    grid on
    legend({'$neutral steering$','$-DeltaAlpha$','$tangent$'});
    hold off;


    % Fitting of the NON LINEAR ZONE
    

figure('Name','Handling diagram NON LINEAR','NumberTitle','off'), clf

    % Cut vectors
    cut_value_start_nl = cut_value_end; %Selection of the starting linearizing point is the end of the previous linear zone (Normalized acceleration value)
    index_start_nl = find((Ay/g) > cut_value_start_nl);
    cut_index_start_nl = index_start_nl(1) - 1; %Selection of the starting linearizing point (index value)

    
    % the ending value correspond directly to the final point of the curvature
    cut_index_end_nl = numel(Ay)-1; %Selection of the ending linearizing point (index value)

  
    fprintf('Il cut index  start vale %d \n', cut_index_start_nl);
    fprintf('Il cut index  end vale %d \n', cut_index_end_nl);
   

    % Fitting to compute the tangent (LINEAR ZONE)
    x_cut_nl = Ay(cut_index_start_nl:cut_index_end_nl)/g;
    y_cut_nl = -delta_alpha_dt((cut_index_start_nl-1):(cut_index_end_nl-1));
    
    %Generation of the new curve
    dim = 3;
    p_nl = polyfit(x_cut_nl, y_cut_nl, dim);
    y_nl = polyval(p_nl, x_cut_nl);
   
    for c=1:(dim+1)
        fprintf('Il coefficiente %d della parte NON lineare vale = %f\n',(dim+1-c), p_nl(c));
    end

    % Creazione della curve
    plot(Ay/g, zeros(size(Ay)), 'Color', color('yellow'),'LineWidth',2);
    hold on;

    plot(Ay/g, -delta_alpha_dt(2:end),'Color',color('blue'),'LineWidth',3);
    plot(Ay/g, y, 'Color',color('green'),'LineWidth',2);
    plot(Ay(1:cut_index_end)/g, y(1:cut_index_end), '--', 'Color',color('red'),'LineWidth',2);
    plot(x_cut_nl, y_nl, '--', 'color',[1 0.5 0] , 'LineWidth',2);
    
    title('Handling diagram')
    ylabel('$-Delta Alpha$ [rad]')
    xlabel('$Ay/g$ [-]')
    grid on
    
    legend({'Neutral steering','-DeltaAlpha','Tangent','Linear fitted','Non linear fitted'}, 'Location', 'southwest');
    hold off;



   












    
    % %% Yaw rate gain - Beta gain
    % % -------------------------------
    % yaw_rate_gain_data = Omega./(delta_D*pi/180);
    % 
    % figure('Name','Yaw rate gain vs u','NumberTitle','off')
    % hold on
    % plot(u*3.6,yaw_rate_gain_data,'LineWidth',2)
    % grid on
    % title('$\Omega/\delta_H$ vs $u$');
    % xlabel('u [km/h]');
    % ylabel('$\Omega/\delta$ [1/s]');
    % hold off
    % 
    % %% Beta gain
    % % -------------------------------
    % beta_gain_data = beta./(delta_D*pi/180);
    % 
    % figure('Name','Beta gain vs u','NumberTitle','off')
    % hold on
    % plot(u*3.6,beta_gain_data,'LineWidth',2)
    % grid on
    % title('$\beta/\delta_H$ vs $u$');
    % xlabel('u [km/h]');
    % ylabel('$\beta/\delta_H$');
    % hold off

end
    
