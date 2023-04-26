%% Assignment 1 - Tyre fitting
% Team 6: Consalvi Natale - Pettene Mattia - Zumerle Matteo

%% Initialization
clc;
close all;
clear;

addpath('dataset/');
addpath('tyre_lib/');

% Tyre geometric data:
diameter = 18*2.56; %
Fz0 = 220;   % [N] nominal load
R0  = diameter/2/100; % [m] get from nominal load R0 (m)

% Constants for angle conversions
to_rad = pi/180;
to_deg = 180/pi;

% Last figure of the pure longitudinal part (to add this new main to the
% other)
last_fig_FX0 = 0;

%% Loading of dataset - first cutting process

data_set_path = 'dataset/';

data_set = 'Hoosier_B1464run23'; % pure lateral forces

fprintf('Loading dataset: ')

switch data_set
  case 'Hoosier_B1464run23'
      fprintf('for pure lateral force analysis.')
      load ([data_set_path, data_set]); % pure lateral

  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion (at the higher pressure)
cut_start_pl = 31350;
cut_end_pl   = 54500;


smpl_range_pl = cut_start_pl:cut_end_pl;

fprintf('\ncompleted!')

%% Plot of the raw data

figure ('Name','Raw dataset', 'NumberTitle',1+ last_fig_FX0)
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start_pl cut_start_pl],y_range,'--r')
plot([cut_end_pl cut_end_pl],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start_pl cut_start_pl],y_range,'--r')
plot([cut_end_pl cut_end_pl],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start_pl cut_start_pl],y_range,'--r')
plot([cut_end_pl cut_end_pl],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start_pl cut_start_pl],y_range,'--r')
plot([cut_end_pl cut_end_pl],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start_pl cut_start_pl],y_range,'--r')
plot([cut_end_pl cut_end_pl],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start_pl cut_start_pl],y_range,'--r')
plot([cut_end_pl cut_end_pl],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')
%% Selection of different ranges and plotting
% for camber angle and vertical load

vec_samples_pl = 1:1:length(smpl_range_pl);

tyre_data_pl = table(); % create empty table
% store raw data in table
tyre_data_pl.SL =  SL(smpl_range_pl);
tyre_data_pl.SA = -SA(smpl_range_pl)*to_rad;    % SAE -> Adapted SAE
tyre_data_pl.FZ = -FZ(smpl_range_pl);           % SAE -> Adapted SAE
tyre_data_pl.FX =  FX(smpl_range_pl);
tyre_data_pl.FY =  FY(smpl_range_pl);   
tyre_data_pl.MZ =  MZ(smpl_range_pl);
tyre_data_pl.IA =  IA(smpl_range_pl)*to_rad;

% Extract points at constant camber angle
GAMMA_tol_pl = 0.05*to_rad;
idx_pl.GAMMA_0 = 0.0*to_rad-GAMMA_tol_pl < tyre_data_pl.IA & tyre_data_pl.IA < 0.0*to_rad+GAMMA_tol_pl;
idx_pl.GAMMA_1 = 1.0*to_rad-GAMMA_tol_pl < tyre_data_pl.IA & tyre_data_pl.IA < 1.0*to_rad+GAMMA_tol_pl;
idx_pl.GAMMA_2 = 2.0*to_rad-GAMMA_tol_pl < tyre_data_pl.IA & tyre_data_pl.IA < 2.0*to_rad+GAMMA_tol_pl;
idx_pl.GAMMA_3 = 3.0*to_rad-GAMMA_tol_pl < tyre_data_pl.IA & tyre_data_pl.IA < 3.0*to_rad+GAMMA_tol_pl;
idx_pl.GAMMA_4 = 4.0*to_rad-GAMMA_tol_pl < tyre_data_pl.IA & tyre_data_pl.IA < 4.0*to_rad+GAMMA_tol_pl;

GAMMA_0_pl  = tyre_data_pl( idx_pl.GAMMA_0, : );
GAMMA_1_pl  = tyre_data_pl( idx_pl.GAMMA_1, : );
GAMMA_2_pl  = tyre_data_pl( idx_pl.GAMMA_2, : );
GAMMA_3_pl  = tyre_data_pl( idx_pl.GAMMA_3, : );
GAMMA_4_pl  = tyre_data_pl( idx_pl.GAMMA_4, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol_pl = 100;
idx_pl.FZ_220  = 220-FZ_tol_pl < tyre_data_pl.FZ & tyre_data_pl.FZ < 220+FZ_tol_pl;
idx_pl.FZ_440  = 440-FZ_tol_pl < tyre_data_pl.FZ & tyre_data_pl.FZ < 440+FZ_tol_pl;
idx_pl.FZ_700  = 700-FZ_tol_pl < tyre_data_pl.FZ & tyre_data_pl.FZ < 700+FZ_tol_pl;
idx_pl.FZ_900  = 900-FZ_tol_pl < tyre_data_pl.FZ & tyre_data_pl.FZ < 900+FZ_tol_pl;
idx_pl.FZ_1120 = 1120-FZ_tol_pl < tyre_data_pl.FZ & tyre_data_pl.FZ < 1120+FZ_tol_pl;
FZ_220_pl  = tyre_data_pl( idx_pl.FZ_220, : );
FZ_440_pl  = tyre_data_pl( idx_pl.FZ_440, : );
FZ_700_pl  = tyre_data_pl( idx_pl.FZ_700, : );
FZ_900_pl  = tyre_data_pl( idx_pl.FZ_900, : );
FZ_1120_pl = tyre_data_pl( idx_pl.FZ_1120, : );

% Plot
figure('Name','Ranges selection', 'NumberTitle', 2 + last_fig_FX0)
tiledlayout(2,1)

ax_list(1) = nexttile;
plot(tyre_data_pl.IA*to_deg)
hold on
plot(vec_samples_pl(idx_pl.GAMMA_0),GAMMA_0_pl.IA*to_deg,'.');
plot(vec_samples_pl(idx_pl.GAMMA_1),GAMMA_1_pl.IA*to_deg,'.');
plot(vec_samples_pl(idx_pl.GAMMA_2),GAMMA_2_pl.IA*to_deg,'.');
plot(vec_samples_pl(idx_pl.GAMMA_3),GAMMA_3_pl.IA*to_deg,'.');
plot(vec_samples_pl(idx_pl.GAMMA_4),GAMMA_4_pl.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data_pl.FZ)
hold on
plot(vec_samples_pl(idx_pl.FZ_220),FZ_220_pl.FZ,'.');
plot(vec_samples_pl(idx_pl.FZ_440),FZ_440_pl.FZ,'.');
plot(vec_samples_pl(idx_pl.FZ_700),FZ_700_pl.FZ,'.');
plot(vec_samples_pl(idx_pl.FZ_900),FZ_900_pl.FZ,'.');
plot(vec_samples_pl(idx_pl.FZ_1120),FZ_1120_pl.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')
hold off

%% Pure conditions range and plotting
% choose the range with: longitudinal slip = 0, camber angle = 0, vertical
% load = Fz0 = 220N (obv within the higher pressure dataset)

[TData0_pl, ~] = intersect_table_data( GAMMA_0_pl, FZ_220_pl );

% % Plot
% figure('Name','Pure conditions range', 'NumberTitle', 3 + last_fig_FX0)
% tiledlayout(2,1)
% 
% ax_list(1) = nexttile;
% plot(TData0_pl.IA*to_deg)
% hold on
% title('Camber angle')
% xlabel('Samples [-]')
% ylabel('[deg]')
% 
% ax_list(2) = nexttile;
% plot(TData0_pl.FZ)
% hold on
% title('Vertical force')
% xlabel('Samples [-]')
% ylabel('[N]')

figure('Name','Pure conditions range', 'NumberTitle', 3 + last_fig_FX0)
plot_selected_data(TData0_pl);

%% Initialization phase for fitting the pure conditions

tyre_coeffs_pl = initialise_tyre_data_ply(R0, Fz0);

%% Pure conditions fitting: Fz0 = 220N, gamma = 0
% Fit the coefficients {pCy1, pDy1, pEy1, pHy1, pKy1, pKy2, pVy1}

%FZ0 = mean(TData0.FZ);

zeros_vec_pl = zeros(size(TData0_pl.SA));
ones_vec_pl  = ones(size(TData0_pl.SA));

FY0_guess = MF96_FY0_vec(zeros_vec_pl, TData0_pl.SA , zeros_vec_pl, tyre_coeffs_pl.FZ0*ones_vec_pl, tyre_coeffs_pl);

% Check lateral pure force guess
figure('Name','FY0 guess', 'NumberTitle', 4 + last_fig_FX0)
plot(TData0_pl.SA,TData0_pl.FY,'.')
hold on
plot(TData0_pl.SA,FY0_guess,'.')
hold off

% Guess values for parameters to be optimised
%       [pCy1 pDy1 pEy1 pHy1  pKy1  pKy2  pVy1]
P0_pl = [ 1, 1,	1,	1,	1,	1,	1 ]; 
%P0_pl = [ 0.643600000000000,	-4.96800000000000,	1.04000000000000,	0.00430000000000000,	-150.090000000000,	-4.48100000000000,	-0.103700000000000 ]; 

% Limits for parameters to be optimised
lb_pl = [ ];
ub_pl = [ ];


ALPHA_vec = TData0_pl.SA;
FY_vec    = TData0_pl.FY;

% check guess
SA_vec = (-12.5*to_rad):0.001:(12.5*to_rad);
FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec , zeros(size(SA_vec)), ...
                              mean(TData0_pl.FZ).*ones(size(SA_vec)),tyre_coeffs_pl);


% Minimize the residual varying X. It is an unconstrained minimization problem 

[P_fz_nom_pl,fval,exitflag] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,mean(TData0_pl.FZ), tyre_coeffs_pl),...
                               P0_pl,[],[],[],[],lb_pl,ub_pl);

% Update tyre data with new optimal values                             
tyre_coeffs_pl.pCy1 = P_fz_nom_pl(1) ; % 1
tyre_coeffs_pl.pDy1 = P_fz_nom_pl(2) ;  
tyre_coeffs_pl.pEy1 = P_fz_nom_pl(3) ;
tyre_coeffs_pl.pHy1 = P_fz_nom_pl(4) ;
tyre_coeffs_pl.pKy1 = P_fz_nom_pl(5) ; 
tyre_coeffs_pl.pKy2 = P_fz_nom_pl(6) ;
tyre_coeffs_pl.pVy1 = P_fz_nom_pl(7) ;

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec , zeros(size(SA_vec)), ...
                              mean(TData0_pl.FZ).*ones(size(SA_vec)),tyre_coeffs_pl);

% Result of the fitting FY0 in the pure conditions
figure('Name','Fy0(Fz0)','NumberTitle', 5 + last_fig_FX0)
plot(TData0_pl.SA,TData0_pl.FY,'o')
hold on
plot(SA_vec,FY0_fz_nom_vec,'.','LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')

%% Fit coefficients with variable load
% extract data with variable load and camber angle equal to 0
TDataDFz_pl = GAMMA_0_pl;

% figure('Name','Variable Fz range', 'NumberTitle', 6 + last_fig_FX0)
% plot_selected_data(TDataDFz);

smpl_range_pl_dFz = size(TDataDFz_pl);
vec_samples_pl_dFz = 1:1:smpl_range_pl_dFz;

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol_pl_dFz = 100;
idx_pl_dFz.FZ_220  = 220-FZ_tol_pl_dFz < TDataDFz_pl.FZ & TDataDFz_pl.FZ < 220+FZ_tol_pl_dFz;
idx_pl_dFz.FZ_440  = 440-FZ_tol_pl_dFz < TDataDFz_pl.FZ & TDataDFz_pl.FZ < 440+FZ_tol_pl_dFz;
idx_pl_dFz.FZ_700  = 700-FZ_tol_pl_dFz < TDataDFz_pl.FZ & TDataDFz_pl.FZ < 700+FZ_tol_pl_dFz;
idx_pl_dFz.FZ_900  = 900-FZ_tol_pl_dFz < TDataDFz_pl.FZ & TDataDFz_pl.FZ < 900+FZ_tol_pl_dFz;
idx_pl_dFz.FZ_1120 = 1120-FZ_tol_pl_dFz < TDataDFz_pl.FZ & TDataDFz_pl.FZ < 1120+FZ_tol_pl_dFz;
FZ_220_pl_dFz  = TDataDFz_pl( idx_pl_dFz.FZ_220, : );
FZ_440_pl_dFz  = TDataDFz_pl( idx_pl_dFz.FZ_440, : );
FZ_700_pl_dFz  = TDataDFz_pl( idx_pl_dFz.FZ_700, : );
FZ_900_pl_dFz  = TDataDFz_pl( idx_pl_dFz.FZ_900, : );
FZ_1120_pl_dFz = TDataDFz_pl( idx_pl_dFz.FZ_1120, : );

% Plot
figure('Name','Considered dataset', 'NumberTitle', 6 + last_fig_FX0)
tiledlayout(2,1)

ax_list_3(1) = nexttile;
plot(TDataDFz_pl.IA*to_deg)
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_3(2) = nexttile;
plot(TDataDFz_pl.FZ)
hold on
plot(vec_samples_pl_dFz(idx_pl_dFz.FZ_220),FZ_220_pl_dFz.FZ,'.');
plot(vec_samples_pl_dFz(idx_pl_dFz.FZ_440),FZ_440_pl_dFz.FZ,'.');
plot(vec_samples_pl_dFz(idx_pl_dFz.FZ_700),FZ_700_pl_dFz.FZ,'.');
plot(vec_samples_pl_dFz(idx_pl_dFz.FZ_900),FZ_900_pl_dFz.FZ,'.');
plot(vec_samples_pl_dFz(idx_pl_dFz.FZ_1120),FZ_1120_pl_dFz.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')
hold off


zeros_vec_pl = zeros(size(TDataDFz_pl.SA));
ones_vec_pl  = ones(size(TDataDFz_pl.SA));

% FY0_guess_dFz = MF96_FY0_vec(zeros_vec_pl , TDataDFz_pl.SA, zeros_vec_pl, tyre_coeffs_pl.FZ0*ones_vec_pl, tyre_coeffs_pl);
% 
% % check guess 
% figure('Name','Check guess for variable Fz with SA DATA', 'NumberTitle', 7 + last_fig_FX0)
% plot(TDataDFz_pl.SA,TDataDFz_pl.FY,'.')
% hold on
% plot(TDataDFz_pl.SA,FY0_guess_dFz,'.')


% Guess values for parameters to be optimised
%    [pDy2 pEy2 pHy2 pVy2] 
P0_pl_dFz =[ 0,0,0,0 ];   

% Limits for parameters to be optimised
%    [pDy2 pEy2 pHy2 pVy2] 
lb_dFz = [ ];
ub_dFz = [ ];


ALPHA_vec_dFz = TDataDFz_pl.SA;
FY_vec_dFz    = TDataDFz_pl.FY;
FZ_vec_dFz    = TDataDFz_pl.FZ;

% check guess

FY0_dfz_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              mean(FZ_220_pl_dFz.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);

figure('Name','Check guess for variable Fz (based on pure conditions)', 'NumberTitle', 8 + last_fig_FX0)
plot(ALPHA_vec_dFz,FY_vec_dFz,'.')
hold on
plot(SA_vec,FY0_dfz_vec,'.')

% Resitual minimization
[P_dfz_pl,fval,exitflag] = fmincon(@(P_pl)resid_pure_Fy_varFz(P_pl,FY_vec_dFz, ALPHA_vec_dFz,0,FZ_vec_dFz, tyre_coeffs_pl),...
                               P0_pl_dFz,[],[],[],[],lb_dFz,ub_dFz);
disp(exitflag)

% Change tyre data with new optimal values                             
tyre_coeffs_pl.pDy2 = P_dfz_pl(1);
tyre_coeffs_pl.pEy2 = P_dfz_pl(2);
tyre_coeffs_pl.pHy2 = P_dfz_pl(3);
tyre_coeffs_pl.pVy2 = P_dfz_pl(4);


res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz_pl , FY_vec_dFz,SA_vec, 0 , FZ_vec_dFz, tyre_coeffs_pl);

tmp_zeros_dFz = zeros(size(SA_vec));
tmp_ones_dFz = ones(size(SA_vec));

FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_220_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_440_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_700_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_900_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_1120_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);


figure('Name','Fy0(Fz)','NumberTitle', 9 + last_fig_FX0)
plot(TDataDFz_pl.SA*to_deg,TDataDFz_pl.FY,'o')
hold on
plot(SA_vec*to_deg,FY0_fz_var_vec1,'.','LineWidth',2)
plot(SA_vec*to_deg,FY0_fz_var_vec2,'.','LineWidth',2)
plot(SA_vec*to_deg,FY0_fz_var_vec3,'.','LineWidth',2)
plot(SA_vec*to_deg,FY0_fz_var_vec4,'.','LineWidth',2)
plot(SA_vec*to_deg,FY0_fz_var_vec5,'.','LineWidth',2)
legend({'data', '$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'}, 'Location','eastoutside');
xlabel('$\alpha$ [-]')
ylabel('$F_{y0}$ [N]')

% %% Fit coefficient with variable camber
% % extract data with variable load
% [TDataGamma_pl, ~] = intersect_table_data( SA_0, FZ_220 );
% 
% % Fit the coeffs { pDx3}
% 
% % Guess values for parameters to be optimised
% P0 = [0]; 
% 
% % NOTE: many local minima => limits on parameters are fundamentals
% % Limits for parameters to be optimised
% % 1< pCx1 < 2 
% % 0< pEx1 < 1 
% %lb = [0, 0,  0, 0,  0,  0,  0];
% %ub = [2, 1e6,1, 1,1e1,1e2,1e2];
% lb = [];
% ub = [];
% 
% 
% zeros_vec = zeros(size(TDataGamma_pl.SL));
% ones_vec  = ones(size(TDataGamma_pl.SL));
% 
% KAPPA_vec = TDataGamma_pl.SL;
% GAMMA_vec = TDataGamma_pl.IA; 
% FX_vec    = TDataGamma_pl.FX;
% FZ_vec    = TDataGamma_pl.FZ;
% 
% figure()
% plot(KAPPA_vec,FX_vec);
% 
% 
% % LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% % is an unconstrained minimization problem 
% [P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
%                                P0,[],[],[],[],lb,ub);
% 
% % Change tyre data with new optimal values                             
% tyre_coeffs.pDx3 = P_varGamma(1) ; % 1
% 
% FX0_varGamma_vec = MF96_FX0_vec(KAPPA_vec,zeros_vec , GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
% 
% figure('Name','Fx0 vs Gamma')
% plot(KAPPA_vec,TDataGamma_pl.FX,'o')
% hold on
% plot(KAPPA_vec,FX0_varGamma_vec,'-')
% xlabel('$\kappa$ [-]')
% ylabel('$F_{x0}$ [N]')
% % Calculate the residuals with the optimal solution found above
% res_Fx0_varGamma  = resid_pure_Fx_varGamma(P_varGamma,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);
% 
% % R-squared is 
% % 1-SSE/SST
% % SSE/SST = res_Fx0_nom
% 
% % SSE is the sum of squared error,  SST is the sum of squared total
% fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);
% 
% 
% [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% % 
% fprintf('Bx      = %6.3f\n',Bx);
% fprintf('Cx      = %6.3f\n',Cx);
% fprintf('mux      = %6.3f\n',Dx/tyre_coeffs.FZ0);
% fprintf('Ex      = %6.3f\n',Ex);
% fprintf('SVx     = %6.3f\n',SVx);
% fprintf('kappa_x = %6.3f\n',kappa__x);
% fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);
% 
% % % Longitudinal stiffness
% % Kx_vec = zeros(size(load_vec));
% % for i = 1:length(load_vec)
% %   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(i), tyre_data);
% %   Kx_vec(i) = Bx*Cx*Dx/tyre_data.Fz0;
% % end
% % 
% % figure('Name','Kx vs Fz')
% % plot(load_vec,Kx_vec,'o-')
