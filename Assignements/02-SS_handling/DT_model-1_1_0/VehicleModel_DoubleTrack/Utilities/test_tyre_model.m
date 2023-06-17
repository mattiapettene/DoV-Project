%% Test Tyre model

load_MF96_tyre_data;

alpha_swipe = deg2rad(-15:0.1:15);
kappa_swipe = -0.3:0.01:0.3;

alpha_vec = deg2rad(linspace(0, 10, 5))';
phi_vec = deg2rad(linspace(0, 5, 5))';
order_magn = floor(log(abs(tyre_data_f.Fz0))./log(10));
fz_vec = linspace(tyre_data_f.Fz0 - 10 ^ order_magn, tyre_data_f.Fz0 + 10 ^ order_magn, 5);
kappa_vec = linspace(0, 0.4, 5)';
fz_nom = tyre_data_f.Fz0;


%% Lateral Force
% Pure Fy var Fz
figure
hold on
tmp_var_vec = fz_vec;
for jj = 1:length(tmp_var_vec)
    plot(alpha_swipe,...
        forloop(alpha_swipe,@(alpha)MF96_FY0(0,alpha,0,tmp_var_vec(jj),tyre_data_f)),...
        "DisplayName",strcat("$F_z$ = ",num2str(tmp_var_vec(jj))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\alpha$ [rad]'); ylabel('$F_y$ [N]')
title('Pure $F_y$ varying $F_z$')

% Pure Fy var Fz
figure
hold on
tmp_var_vec = phi_vec;
for jj = 1:length(tmp_var_vec)
    plot(alpha_swipe,...
        forloop(alpha_swipe,@(alpha)MF96_FY0(0,alpha,tmp_var_vec(jj),fz_nom,tyre_data_f)),...
        "DisplayName",strcat("$\gamma$ = ",num2str(rad2deg(tmp_var_vec(jj)))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\alpha$ [rad]'); ylabel('$F_y$ [N]')
title('Pure $F_y$ varying $\gamma$')

% Fy var kappa
figure
hold on
tmp_var_vec = kappa_vec;
for jj = 1:length(tmp_var_vec)
    plot(alpha_swipe,...
        forloop(alpha_swipe,@(alpha)MF96_FY(tmp_var_vec(jj),alpha,0,fz_nom,tyre_data_f)),...
        "DisplayName",strcat("$\kappa$ = ",num2str(rad2deg(tmp_var_vec(jj)))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\alpha$ [rad]'); ylabel('$F_y$ [N]')
title('Combined $F_y$ varying $\kappa$')

%% Longitudinal Force
% Pure Fx var Fz
figure
hold on
tmp_var_vec = fz_vec;
for jj = 1:length(tmp_var_vec)
    plot(kappa_swipe,...
        forloop(kappa_swipe,@(kappa)MF96_FX0(kappa,0,0,tmp_var_vec(jj),tyre_data_f)),...
        "DisplayName",strcat("$F_z$ = ",num2str(tmp_var_vec(jj))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\kappa$ [-]'); ylabel('$F_x$ [N]')
title('Pure $F_x$ varying $F_z$')

% Pure Fy var Fz
figure
hold on
tmp_var_vec = phi_vec;
for jj = 1:length(tmp_var_vec)
    plot(kappa_swipe,...
        forloop(kappa_swipe,@(kappa)MF96_FX0(kappa,0,tmp_var_vec(jj),fz_nom,tyre_data_f)),...
        "DisplayName",strcat("$\gamma$ = ",num2str(rad2deg(tmp_var_vec(jj)))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\kappa$ [-]'); ylabel('$F_x$ [N]')
title('Pure $F_x$ varying $\gamma$')

% Fy var kappa
figure
hold on
tmp_var_vec = alpha_vec;
for jj = 1:length(tmp_var_vec)
    plot(kappa_swipe,...
        forloop(kappa_swipe,@(kappa)MF96_FX(kappa,tmp_var_vec(jj),0,fz_nom,tyre_data_f)),...
        "DisplayName",strcat("$\alpha$ = ",num2str(rad2deg(tmp_var_vec(jj)))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\kappa$ [-]'); ylabel('$F_x$ [N]')
title('Combined $F_x$ varying $\alpha$')

%% Aligning torque

%% Lateral Force
% Pure Fy var Fz
figure
hold on
tmp_var_vec = fz_vec;
for jj = 1:length(tmp_var_vec)
    plot(alpha_swipe,...
        forloop(alpha_swipe,@(alpha)MF96_MZ0(0,alpha,0,tmp_var_vec(jj),tyre_data_f)),...
        "DisplayName",strcat("$F_z$ = ",num2str(tmp_var_vec(jj))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\alpha$ [rad]'); ylabel('$F_y$ [N]')
title('Pure $F_y$ varying $F_z$')

% Pure Fy var Fz
figure
hold on
tmp_var_vec = phi_vec;
for jj = 1:length(tmp_var_vec)
    plot(alpha_swipe,...
        forloop(alpha_swipe,@(alpha)MF96_MZ0(0,alpha,tmp_var_vec(jj),fz_nom,tyre_data_f)),...
        "DisplayName",strcat("$\gamma$ = ",num2str(rad2deg(tmp_var_vec(jj)))),...
        "LineWidth",2);
end
hold off
legend
grid
xlabel('$\alpha$ [rad]'); ylabel('$F_y$ [N]')
title('Pure $F_y$ varying $\gamma$')







%% Aux functions

function res = forloop(loop_vec,fun_handle)
for k = 1:length(loop_vec)
    res(k,1) = fun_handle(loop_vec(k));
end
end