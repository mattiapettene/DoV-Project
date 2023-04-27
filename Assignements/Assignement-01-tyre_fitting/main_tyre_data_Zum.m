%% Assignment 1 - Tyre fitting
% Team 6: Consalvi Natale - Pettene Mattia - Zumerle Matteo

%% Initialization
clc;
close all;
clear;

set(0,'DefaultFigureWindowStyle','docked');

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
tiledlayout(7,1)

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

ax_list(7) = nexttile; y_range = [min(min(MZ),0) round(max(MZ)*1.1)];
plot(MZ)
hold on
plot([cut_start_pl cut_start_pl],y_range,'--r')
plot([cut_end_pl cut_end_pl],y_range,'--r')
title('Aligning moment Mz')
xlabel('Samples [-]')
ylabel('[Nm]')

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

ax_list_2(1) = nexttile;
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

ax_list_2(2) = nexttile;
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
linkaxes(ax_list_2,'x')
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
linkaxes(ax_list_3,'x')

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

%% Fit coefficient with variable camber
% extract data with variable load
TDataGamma_pl = FZ_220_pl;

smpl_range_pl_dgamma = size(TDataGamma_pl);
vec_samples_pl_dgamma = 1:1:smpl_range_pl_dgamma;

% Extract points at constant camber and plot
GAMMA_tol_pl_dgamma = 0.05*to_rad;
idx_pl_dgamma.GAMMA_0 = 0.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 0.0*to_rad+GAMMA_tol_pl_dgamma;
idx_pl_dgamma.GAMMA_1 = 1.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 1.0*to_rad+GAMMA_tol_pl_dgamma;
idx_pl_dgamma.GAMMA_2 = 2.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 2.0*to_rad+GAMMA_tol_pl_dgamma;
idx_pl_dgamma.GAMMA_3 = 3.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 3.0*to_rad+GAMMA_tol_pl_dgamma;
idx_pl_dgamma.GAMMA_4 = 4.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 4.0*to_rad+GAMMA_tol_pl_dgamma;

GAMMA_0_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_0, : );
GAMMA_1_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_1, : );
GAMMA_2_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_2, : );
GAMMA_3_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_3, : );
GAMMA_4_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_4, : );

% Plot
figure('Name','Considered dataset for variable camber', 'NumberTitle', 10 + last_fig_FX0)
tiledlayout(3,1)
ax_list_4(1) = nexttile;
plot(TDataGamma_pl.IA*to_deg)
hold on
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_0),GAMMA_0_dgamma.IA*to_deg,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_1),GAMMA_1_dgamma.IA*to_deg,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_2),GAMMA_2_dgamma.IA*to_deg,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_3),GAMMA_3_dgamma.IA*to_deg,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_4),GAMMA_4_dgamma.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')
hold off

ax_list_4(2) = nexttile;
plot(TDataGamma_pl.FY)
hold on
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_0),GAMMA_0_dgamma.FY,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_1),GAMMA_1_dgamma.FY,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_2),GAMMA_2_dgamma.FY,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_3),GAMMA_3_dgamma.FY,'.');
plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_4),GAMMA_4_dgamma.FY,'.');
title('Lateral force')
xlabel('Samples [-]')
ylabel('[N]')
hold off

ax_list_4(3) = nexttile;
plot(TDataGamma_pl.FZ)
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')
linkaxes(ax_list_4,'x')

% Fit the coeffs {pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4}

% Guess values for parameters to be optimised
% Fit the coeffs { pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4 }
% Guess values for parameters to be optimised
P0_pl_dgamma = [0.51e1, 1.88, -0.42e1, -2.04, 0.13e1, -0.29e1, -0.28e1];

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%lb = [0, 0,  0, 0,  0,  0,  0];
%ub = [2, 1e6,1, 1,1e1,1e2,1e2];
lb_dgamma = [5,1,-100,-2.5,0,-100,-100];
ub_dgamma = [100,100,0,-0.1,100,0,0];

zeros_vec_dgamma = zeros(size(TDataGamma_pl.IA));
ones_vec_dgamma  = ones(size(TDataGamma_pl.IA));

ALPHA_vec_dgamma = TDataGamma_pl.SA; 
GAMMA_vec_dgamma = TDataGamma_pl.IA; 
FY_vec_dgamma    = TDataGamma_pl.FY;
FZ_vec_dgamma    = TDataGamma_pl.FZ;

figure('Name','Non so cosa sia', 'NumberTitle', 11 + last_fig_FX0)
plot(ALPHA_vec_dgamma,FY_vec_dgamma);


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec_dgamma, ALPHA_vec_dgamma,GAMMA_vec_dgamma,tyre_coeffs_pl.FZ0, tyre_coeffs_pl),...
                               P0_pl_dgamma,[],[],[],[],lb_dgamma,ub_dgamma);

% Change tyre data with new optimal values                             
tyre_coeffs_pl.pDy3 = P_varGamma(1);  
tyre_coeffs_pl.pEy3 = P_varGamma(2); 
tyre_coeffs_pl.pEy4 = P_varGamma(3); 
tyre_coeffs_pl.pHy3 = P_varGamma(4); 
tyre_coeffs_pl.pKy3 = P_varGamma(5); 
tyre_coeffs_pl.pVy3 = P_varGamma(6); 
tyre_coeffs_pl.pVy4 = P_varGamma(7); 

FY0_varGamma_vec = MF96_FY0_vec(zeros_vec_dgamma, ALPHA_vec_dgamma , GAMMA_vec_dgamma, tyre_coeffs_pl.FZ0*ones_vec_dgamma,tyre_coeffs_pl);

figure('Name','Fx0 vs Gamma', 'NumberTitle', 12 + last_fig_FX0)
plot(ALPHA_vec_dgamma,TDataGamma_pl.FY,'o')
hold on
plot(ALPHA_vec_dgamma,FY0_varGamma_vec,'-')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

tmp_zeros_dgamma = zeros(size(SA_vec));
tmp_ones_dgamma = ones(size(SA_vec));

FY0_gamma_var_vec1 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_0_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
FY0_gamma_var_vec2 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_1_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
FY0_gamma_var_vec3 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_2_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
FY0_gamma_var_vec4 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_3_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
FY0_gamma_var_vec5 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_4_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);


figure('Name','Last fig','NumberTitle', 13 + last_fig_FX0)
hold on
plot(GAMMA_0_dgamma.SA*to_deg,GAMMA_0_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','y',
plot(GAMMA_1_dgamma.SA*to_deg,GAMMA_1_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','c',
plot(GAMMA_2_dgamma.SA*to_deg,GAMMA_2_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','m',
plot(GAMMA_3_dgamma.SA*to_deg,GAMMA_3_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','b',
plot(GAMMA_4_dgamma.SA*to_deg,GAMMA_4_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','r',
plot(SA_vec*to_deg,FY0_gamma_var_vec1,'-s','LineWidth',2,'MarkerSize',1)
plot(SA_vec*to_deg,FY0_gamma_var_vec2,'-s','LineWidth',2,'MarkerSize',1)
plot(SA_vec*to_deg,FY0_gamma_var_vec3,'-s','LineWidth',2,'MarkerSize',1)
plot(SA_vec*to_deg,FY0_gamma_var_vec4,'-s','LineWidth',2,'MarkerSize',1)
plot(SA_vec*to_deg,FY0_gamma_var_vec5,'-s','LineWidth',2,'MarkerSize',1)
legend({'$ \gamma_0 = 0 deg $','$ \gamma_1 = 1 deg $','$ \gamma_2 = 2 deg $','$ \gamma_3 = 3 deg$','$ \gamma_4 = 4 deg $', 'Fy($\gamma_0$)','Fy($\gamma_1$)','Fy($\gamma_2$)','Fy($\gamma_3$)','Fy($\gamma_4$)'}, 'Location','eastoutside');
xlabel('$\alpha$ [-]')
ylabel('$F_{y0}$ [N]')
legend({'Fy($\gamma_0$)','Fy($\gamma_1$)','Fy($\gamma_2$)','Fy($\gamma_3$)','Fy($\gamma_4$)'}, 'Location','eastoutside');


% Calculate the residuals with the optimal solution found above
res_Fy0_dgamma  = resid_pure_Fy_varGamma(P_varGamma,FY_vec_dgamma, ALPHA_vec_dgamma,GAMMA_vec_dgamma,tyre_coeffs_pl.FZ0, tyre_coeffs_pl);

% % R-squared is 
% % 1-SSE/SST
% % SSE/SST = res_Fx0_nom
% 
% % SSE is the sum of squared error,  SST is the sum of squared total
% fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);
% 
% 
% [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec_dgamma(3), tyre_coeffs.FZ0, tyre_coeffs);
% % 
% fprintf('Bx      = %6.3f\n',Bx);
% fprintf('Cx      = %6.3f\n',Cx);
% fprintf('mux      = %6.3f\n',Dx/tyre_coeffs.FZ0);
% fprintf('Ex      = %6.3f\n',Ex);
% fprintf('SVx     = %6.3f\n',SVx);
% fprintf('kappa_x = %6.3f\n',kappa__x);
% fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);

% % Longitudinal stiffness
% Kx_vec = zeros(size(load_vec));
% for i = 1:length(load_vec)
%   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(i), tyre_data);
%   Kx_vec(i) = Bx*Cx*Dx/tyre_data.Fz0;
% end
% 
% figure('Name','Kx vs Fz')
% plot(load_vec,Kx_vec,'o-')

%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------------------------------
%% %% Pure conditions range and plotting

%%Intersect table data
GAMMA_0_mz=GAMMA_0_pl;
[TData0_mz, ~] = intersect_table_data(GAMMA_0_mz, FZ_220_pl);
%% Initialization phase for fitting the pure conditions

tyre_coeffs_mz = initialise_tyre_data_ply(R0, Fz0);

%% Pure conditions fitting: Fz0 = 220N, gamma = 0
%FZ0 = mean(TData0.FZ);

zeros_vec_mz = zeros(size(TData0_mz.MZ));
ones_vec_mz  = ones(size(TData0_mz.MZ));

MZ0_guess = MF96_MZ0_vec(zeros_vec_mz, TData0_mz.SA , zeros_vec_mz, tyre_coeffs_mz.FZ0*ones_vec_mz, tyre_coeffs_mz);

% Check aligning moment guess
figure('Name','MZ0 guess', 'NumberTitle', 4 + last_fig_FX0)
plot(TData0_mz.SA,TData0_mz.MZ,'.')
hold on
plot(TData0_mz.SA,MZ0_guess,'.')
hold off
% %% 
% 
% Guess values for parameters to be optimised
% %      {qBz1, qBz9, qBz10, qCz1, qDz1, qDz6, qEz1, qEz4, qHz1}
%P0_mz = [ 3.8, -0.5, -18, 2, -1.3, -0.01, 0.2, 2.8, 0.5]; 
P0_mz = [1, 1, 1, 1, 1, 1, 1, 1, 1];

% % Limits for parameters to be optimised
lb_mz = [];
ub_mz = [];
% 
% 
ALPHA_vec = TData0_mz.SA;
MZ_vec    = TData0_mz.MZ;
% 
% check guess
SA_vec = (-12.5*to_rad):0.001:(12.5*to_rad);
MZ0_fz_nom_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec , zeros(size(SA_vec)), ...
                               mean(TData0_mz.FZ).*ones(size(SA_vec)),tyre_coeffs_mz);
% 
% 
% % Minimize the residual. It is an unconstrained minimization problem 
% 
[Pmz_fz_nom_mz,fval,exitflag] = fmincon(@(P)resid_pure_Mz(P,MZ_vec, ALPHA_vec,0,mean(TData0_mz.FZ), tyre_coeffs_mz),...
                                P0_mz,[],[],[],[],lb_mz,ub_mz);
display(Pmz_fz_nom_mz)

% Update tyre data with new optimal values                             
tyre_coeffs_mz.qBz1 = Pmz_fz_nom_mz(1) ; % 1
tyre_coeffs_mz.qBz9 = Pmz_fz_nom_mz(2) ;
tyre_coeffs_mz.qBz10 = Pmz_fz_nom_mz(3) ;
tyre_coeffs_mz.qCz1 = Pmz_fz_nom_mz(4) ;
tyre_coeffs_mz.qDz1 = Pmz_fz_nom_mz(5) ;
tyre_coeffs_mz.qDz6 = Pmz_fz_nom_mz(6) ;
tyre_coeffs_mz.qEz1 = Pmz_fz_nom_mz(7) ;
tyre_coeffs_mz.qEz4 = Pmz_fz_nom_mz(8) ;
tyre_coeffs_mz.qHz1 = Pmz_fz_nom_mz(9) ;



MZ0_fz_nom_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec , zeros(size(SA_vec)), ...
                               mean(TData0_mz.FZ).*ones(size(SA_vec)),tyre_coeffs_mz);

% Result of the fitting MZ0 in the pure conditions
figure('Name','Mz0(Fz0)','NumberTitle', 5 + last_fig_FX0)
plot(TData0_mz.SA,TData0_mz.MZ,'o')
hold on
plot(SA_vec,MZ0_fz_nom_vec,'.','LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$M_{y0}$ [N]')
%% Fit coefficients with variable load
%extract data with variable load and camber angle equal to 0
TDataDFz_mz = GAMMA_0_mz;

% figure('Name','Variable Fz range', 'NumberTitle', 6 + last_fig_FX0)
% plot_selected_data(TDataDFz);

smpl_range_mz_dFz = size(TDataDFz_mz);
vec_samples_mz_dFz = 1:1:smpl_range_mz_dFz;

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol_mz_dFz = 100;
idx_mz_dFz.FZ_220  = 220-FZ_tol_mz_dFz < TDataDFz_mz.FZ & TDataDFz_mz.FZ < 220+FZ_tol_mz_dFz;
idx_mz_dFz.FZ_440  = 440-FZ_tol_mz_dFz < TDataDFz_mz.FZ & TDataDFz_mz.FZ < 440+FZ_tol_mz_dFz;
idx_mz_dFz.FZ_700  = 700-FZ_tol_mz_dFz < TDataDFz_mz.FZ & TDataDFz_mz.FZ < 700+FZ_tol_mz_dFz;
idx_mz_dFz.FZ_900  = 900-FZ_tol_mz_dFz < TDataDFz_mz.FZ & TDataDFz_mz.FZ < 900+FZ_tol_mz_dFz;
idx_mz_dFz.FZ_1120 = 1120-FZ_tol_mz_dFz < TDataDFz_mz.FZ & TDataDFz_mz.FZ < 1120+FZ_tol_mz_dFz;
FZ_220_mz_dFz  = TDataDFz_mz( idx_mz_dFz.FZ_220, : );
FZ_440_mz_dFz  = TDataDFz_mz( idx_mz_dFz.FZ_440, : );
FZ_700_mz_dFz  = TDataDFz_mz( idx_mz_dFz.FZ_700, : );
FZ_900_mz_dFz  = TDataDFz_mz( idx_mz_dFz.FZ_900, : );
FZ_1120_mz_dFz = TDataDFz_mz( idx_mz_dFz.FZ_1120, :);

% Plot
figure('Name','Considered dataset', 'NumberTitle', 6 + last_fig_FX0)
tiledlayout(2,1)

ax_list_3(1) = nexttile;
plot(TDataDFz_mz.IA*to_deg)
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_3(2) = nexttile;
plot(TDataDFz_mz.MZ)
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
linkaxes(ax_list_3,'x')

zeros_vec_mz = zeros(size(TDataDFz_mz.SA));
ones_vec_mz  = ones(size(TDataDFz_mz.SA));

MZ0_guess_dFz = MF96_MZ0_vec(zeros_vec_mz , TDataDFz_mz.SA, zeros_vec_mz, tyre_coeffs_mz.FZ0*ones_vec_mz, tyre_coeffs_mz);

% check guess 
figure('Name','Check guess for variable Fz with SA DATA', 'NumberTitle', 7 + last_fig_FX0)
plot(TDataDFz_mz.SA,TDataDFz_mz.MZ,'.')
hold on
plot(TDataDFz_mz.SA,MZ0_guess_dFz,'.')


% Guess values for parameters to be optimised
% [qBz2, qBz3, qDz2, qDz7, qEz2, qEz3, qHz2] 

qBz1  = 7.35347522395674,
qBz2  = -4.71595003290210,
qBz3  = -6.17019496981862,
qBz4  = -0.281057344706274,
qBz5  = -0.281558360232466,
qBz9  = 0,
qBz10 = 0,
qCz1  = 1.67097884314219,
qDz1  = 0.197223329858040,
qDz2  = -0.0549111816154103,
qDz3  = -3.31843325657878,
qDz4  = 42.5363302221465,
qDz6  = 0.00307532725022005,
qDz7  = 0.0171847774709159,
qDz8  = -1.82587609907207,
qDz9  = -0.796929452914326,
qEz1  = 0.515207229031359,
qEz2  = -1.69617162475273,
qEz3  = -3.39063846665826,
qEz4  = 0.164755141111871,
qEz5  = -9.68793770573582,
qHz1  = -0.0106709285011952,
qHz2  = -0.00767395160043653,
qHz3  = 0.720767969228803,
qHz4  = 0.0934620620316319,
R0    = 0.228600000000000,

LT  = 1,
LMR = 1


P0_pl_dFz =[ -4.7, -6.1,  -0.05, 0.01, -1.6, -3.39, -0.007 ];   
% 
% % Limits for parameters to be optimised
lb_mz_dFz = [ ];
ub_mz_dFz = [ ];


ALPHA_vec_dFz = TDataDFz_mz.SA;
MZ_vec_dFz    = TDataDFz_mz.MZ;
FZ_vec_dFz    = TDataDFz_mz.FZ;

% check guess

MZ0_dfz_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              mean(FZ_220_pl_dFz.FZ)*ones(size(SA_vec)),tyre_coeffs_mz);

figure('Name','Check guess for variable Fz (based on pure conditions)', 'NumberTitle', 8 + last_fig_FX0)
plot(ALPHA_vec_dFz,MZ_vec_dFz,'.')
hold on
plot(SA_vec,MZ0_dfz_vec,'.')

% Resitual minimization
[P_dfz_mz,fval,exitflag] = fmincon(@(P_mz)resid_pure_Fy_varFz(P_mz,MZ_vec_dFz, ALPHA_vec_dFz,0,FZ_vec_dFz, tyre_coeffs_mz),...
                               Pmz_fz_nom_mz,[],[],[],[],lb_mz_dFz,ub_mz_dFz);
disp(exitflag)

% Change tyre data with new optimal values  
%[qBz2, qBz3, qDz2, qDz7, qEz2, qEz3, qHz2]
tyre_coeffs_mz.qBz2 = P_dfz_mz(1);
tyre_coeffs_mz.qBz3 = P_dfz_mz(2);
tyre_coeffs_mz.qDz2 = P_dfz_mz(3);
tyre_coeffs_mz.qDz7 = P_dfz_mz(4);
tyre_coeffs_mz.qEz2 = P_dfz_mz(5);
tyre_coeffs_mz.qEz3 = P_dfz_mz(6);
tyre_coeffs_mz.qHz2 = P_dfz_mz(7);


res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz_mz , MZ_vec_dFz,SA_vec, 0 , FZ_vec_dFz, tyre_coeffs_mz);

tmp_zeros_dFz = zeros(size(SA_vec));
tmp_ones_dFz = ones(size(SA_vec));

MZ0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_220_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_mz);
MZ0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_440_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_mz);
MZ0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_700_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_mz);
MZ0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_900_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_mz);
MZ0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_1120_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_mz);


figure('Name','Mz0(Fz)','NumberTitle', 9 + last_fig_FX0)
plot(TDataDFz_mz.SA*to_deg,TDataDFz_mz.MZ,'o')
hold on
plot(SA_vec*to_deg,MZ0_fz_var_vec1,'.','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec2,'.','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec3,'.','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec4,'.','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec5,'.','LineWidth',2)
legend({'data', '$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'}, 'Location','eastoutside');
xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [N]')

% %% Fit coefficient with variable camber
% % extract data with variable load
% TDataGamma_pl = FZ_220_pl;
% 
% smpl_range_pl_dgamma = size(TDataGamma_pl);
% vec_samples_pl_dgamma = 1:1:smpl_range_pl_dgamma;
% 
% % Extract points at constant camber and plot
% GAMMA_tol_pl_dgamma = 0.05*to_rad;
% idx_pl_dgamma.GAMMA_0 = 0.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 0.0*to_rad+GAMMA_tol_pl_dgamma;
% idx_pl_dgamma.GAMMA_1 = 1.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 1.0*to_rad+GAMMA_tol_pl_dgamma;
% idx_pl_dgamma.GAMMA_2 = 2.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 2.0*to_rad+GAMMA_tol_pl_dgamma;
% idx_pl_dgamma.GAMMA_3 = 3.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 3.0*to_rad+GAMMA_tol_pl_dgamma;
% idx_pl_dgamma.GAMMA_4 = 4.0*to_rad-GAMMA_tol_pl_dgamma < TDataGamma_pl.IA & TDataGamma_pl.IA < 4.0*to_rad+GAMMA_tol_pl_dgamma;
% 
% GAMMA_0_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_0, : );
% GAMMA_1_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_1, : );
% GAMMA_2_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_2, : );
% GAMMA_3_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_3, : );
% GAMMA_4_dgamma  = TDataGamma_pl( idx_pl_dgamma.GAMMA_4, : );
% 
% % Plot
% figure('Name','Considered dataset for variable camber', 'NumberTitle', 10 + last_fig_FX0)
% tiledlayout(3,1)
% ax_list_4(1) = nexttile;
% plot(TDataGamma_pl.IA*to_deg)
% hold on
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_0),GAMMA_0_dgamma.IA*to_deg,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_1),GAMMA_1_dgamma.IA*to_deg,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_2),GAMMA_2_dgamma.IA*to_deg,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_3),GAMMA_3_dgamma.IA*to_deg,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_4),GAMMA_4_dgamma.IA*to_deg,'.');
% title('Camber angle')
% xlabel('Samples [-]')
% ylabel('[deg]')
% hold off
% 
% ax_list_4(2) = nexttile;
% plot(TDataGamma_pl.FY)
% hold on
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_0),GAMMA_0_dgamma.FY,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_1),GAMMA_1_dgamma.FY,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_2),GAMMA_2_dgamma.FY,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_3),GAMMA_3_dgamma.FY,'.');
% plot(vec_samples_pl_dgamma(idx_pl_dgamma.GAMMA_4),GAMMA_4_dgamma.FY,'.');
% title('Lateral force')
% xlabel('Samples [-]')
% ylabel('[N]')
% hold off
% 
% ax_list_4(3) = nexttile;
% plot(TDataGamma_pl.FZ)
% title('Vertical force')
% xlabel('Samples [-]')
% ylabel('[N]')
% linkaxes(ax_list_4,'x')
% 
% % Fit the coeffs {pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4}
% 
% % Guess values for parameters to be optimised
% % Fit the coeffs { pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4 }
% % Guess values for parameters to be optimised
% P0_pl_dgamma = [0.51e1, 1.88, -0.42e1, -2.04, 0.13e1, -0.29e1, -0.28e1];
% 
% % NOTE: many local minima => limits on parameters are fundamentals
% % Limits for parameters to be optimised
% % 1< pCx1 < 2 
% % 0< pEx1 < 1 
% %lb = [0, 0,  0, 0,  0,  0,  0];
% %ub = [2, 1e6,1, 1,1e1,1e2,1e2];
% lb_dgamma = [5,1,-100,-2.5,0,-100,-100];
% ub_dgamma = [100,100,0,-0.1,100,0,0];
% 
% zeros_vec_dgamma = zeros(size(TDataGamma_pl.IA));
% ones_vec_dgamma  = ones(size(TDataGamma_pl.IA));
% 
% ALPHA_vec_dgamma = TDataGamma_pl.SA; 
% GAMMA_vec_dgamma = TDataGamma_pl.IA; 
% FY_vec_dgamma    = TDataGamma_pl.FY;
% FZ_vec_dgamma    = TDataGamma_pl.FZ;
% 
% figure('Name','Non so cosa sia', 'NumberTitle', 11 + last_fig_FX0)
% plot(ALPHA_vec_dgamma,FY_vec_dgamma);
% 
% 
% % LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% % is an unconstrained minimization problem 
% [P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec_dgamma, ALPHA_vec_dgamma,GAMMA_vec_dgamma,tyre_coeffs_pl.FZ0, tyre_coeffs_pl),...
%                                P0_pl_dgamma,[],[],[],[],lb_dgamma,ub_dgamma);
% 
% % Change tyre data with new optimal values                             
% tyre_coeffs_pl.pDy3 = P_varGamma(1);  
% tyre_coeffs_pl.pEy3 = P_varGamma(2); 
% tyre_coeffs_pl.pEy4 = P_varGamma(3); 
% tyre_coeffs_pl.pHy3 = P_varGamma(4); 
% tyre_coeffs_pl.pKy3 = P_varGamma(5); 
% tyre_coeffs_pl.pVy3 = P_varGamma(6); 
% tyre_coeffs_pl.pVy4 = P_varGamma(7); 
% 
% FY0_varGamma_vec = MF96_FY0_vec(zeros_vec_dgamma, ALPHA_vec_dgamma , GAMMA_vec_dgamma, tyre_coeffs_pl.FZ0*ones_vec_dgamma,tyre_coeffs_pl);
% 
% figure('Name','Fx0 vs Gamma', 'NumberTitle', 12 + last_fig_FX0)
% plot(ALPHA_vec_dgamma,TDataGamma_pl.FY,'o')
% hold on
% plot(ALPHA_vec_dgamma,FY0_varGamma_vec,'-')
% xlabel('$\kappa$ [-]')
% ylabel('$F_{x0}$ [N]')
% 
% tmp_zeros_dgamma = zeros(size(SA_vec));
% tmp_ones_dgamma = ones(size(SA_vec));
% 
% FY0_gamma_var_vec1 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_0_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
% FY0_gamma_var_vec2 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_1_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
% FY0_gamma_var_vec3 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_2_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
% FY0_gamma_var_vec4 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_3_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
% FY0_gamma_var_vec5 = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_4_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
% 
% 
% figure('Name','Last fig','NumberTitle', 13 + last_fig_FX0)
% hold on
% plot(GAMMA_0_dgamma.SA*to_deg,GAMMA_0_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','y',
% plot(GAMMA_1_dgamma.SA*to_deg,GAMMA_1_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','c',
% plot(GAMMA_2_dgamma.SA*to_deg,GAMMA_2_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','m',
% plot(GAMMA_3_dgamma.SA*to_deg,GAMMA_3_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','b',
% plot(GAMMA_4_dgamma.SA*to_deg,GAMMA_4_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','r',
% plot(SA_vec*to_deg,FY0_gamma_var_vec1,'-s','LineWidth',2,'MarkerSize',1)
% plot(SA_vec*to_deg,FY0_gamma_var_vec2,'-s','LineWidth',2,'MarkerSize',1)
% plot(SA_vec*to_deg,FY0_gamma_var_vec3,'-s','LineWidth',2,'MarkerSize',1)
% plot(SA_vec*to_deg,FY0_gamma_var_vec4,'-s','LineWidth',2,'MarkerSize',1)
% plot(SA_vec*to_deg,FY0_gamma_var_vec5,'-s','LineWidth',2,'MarkerSize',1)
% legend({'$ \gamma_0 = 0 deg $','$ \gamma_1 = 1 deg $','$ \gamma_2 = 2 deg $','$ \gamma_3 = 3 deg$','$ \gamma_4 = 4 deg $', 'Fy($\gamma_0$)','Fy($\gamma_1$)','Fy($\gamma_2$)','Fy($\gamma_3$)','Fy($\gamma_4$)'}, 'Location','eastoutside');
% xlabel('$\alpha$ [-]')
% ylabel('$F_{y0}$ [N]')
% legend({'Fy($\gamma_0$)','Fy($\gamma_1$)','Fy($\gamma_2$)','Fy($\gamma_3$)','Fy($\gamma_4$)'}, 'Location','eastoutside');
% 
% 
% % Calculate the residuals with the optimal solution found above
% res_Fy0_dgamma  = resid_pure_Fy_varGamma(P_varGamma,FY_vec_dgamma, ALPHA_vec_dgamma,GAMMA_vec_dgamma,tyre_coeffs_pl.FZ0, tyre_coeffs_pl);
% 
% % % R-squared is 
% % % 1-SSE/SST
% % % SSE/SST = res_Fx0_nom
% % 
% % % SSE is the sum of squared error,  SST is the sum of squared total
% % fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);
% % 
% % 
% % [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec_dgamma(3), tyre_coeffs.FZ0, tyre_coeffs);
% % % 
% % fprintf('Bx      = %6.3f\n',Bx);
% % fprintf('Cx      = %6.3f\n',Cx);
% % fprintf('mux      = %6.3f\n',Dx/tyre_coeffs.FZ0);
% % fprintf('Ex      = %6.3f\n',Ex);
% % fprintf('SVx     = %6.3f\n',SVx);
% % fprintf('kappa_x = %6.3f\n',kappa__x);
% % fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);
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
% 
% 
