%% -Assignment 1 - Tyre fitting-
% Team 6: Consalvi Natale - Pettene Mattia - Zumerle Matteo

%% --Initialization
% define geometric data of tyre, import the path

clc;
close all;
clear;

set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)

addpath('dataset/');
addpath('tyre_lib/');

% Tyre geometric data:
diameter = 18*2.56; %
Fz0 = 220;   % [N] nominal load
R0  = diameter/2/100; % [m] get from nominal load R0 (m)

% Constants for angle conversions
to_rad = pi/180;
to_deg = 180/pi;

data_set_path = 'dataset/';

%% --Initialization phase for tyre coefficients
tyre_coeffs_pl = initialise_tyre_data_ply(R0, Fz0);


%% --Pure longitudinal force FX0: dataset import

data_set = 'Hoosier_B1464run30'; % pure lateral forces

fprintf('Loading dataset: ')

switch data_set
  case 'Hoosier_B1464run30'
      fprintf('for pure longitudinal force analysis.')
      load ([data_set_path, data_set]); % pure lateral

  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion (at the higher pressure)
  cut_start = 19028;
  cut_end   = 37643;


smpl_range = cut_start:cut_end;

fprintf('\ncompleted!')


%% ---Dataset for pure longitudinal: plot

figure ('Name','FX0: entire raw dataset', 'NumberTitle', 1)
tiledlayout(6,1)

ax_list_x(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list_x(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_x(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_x(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list_x(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list_x(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list_x,'x')

%% ---Higher pressure dataset for pure longitudunal force: table selection and plot
% consider the high pressure region of original dataset (more stable one)

vec_samples = 1:1:length(smpl_range);
tyre_data = table();
% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ = -FZ(smpl_range);
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  FY(smpl_range);
tyre_data.MZ =  MZ(smpl_range);
tyre_data.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );

% The slip angle is varied continuously between -4 and +12° and then
% between -12° and +4° for the pure slip case

% The slip angle is varied step wise for longitudinal slip tests
% 0° , - 3° , -6 °
SA_tol = 0.5*to_rad;
idx.SA_0    =  0-SA_tol          < tyre_data.SA & tyre_data.SA < 0+SA_tol;
idx.SA_3neg = -(3*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -3*to_rad+SA_tol;
idx.SA_6neg = -(6*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -6*to_rad+SA_tol;
SA_0     = tyre_data( idx.SA_0, : );
SA_3neg  = tyre_data( idx.SA_3neg, : );
SA_6neg  = tyre_data( idx.SA_6neg, : );

figure('Name','FX0: higher pressure dataset with regions', 'NumberTitle', 2)
tiledlayout(3,1)

ax_list_2_x(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_2_x(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list_2_x(3) = nexttile;
plot(tyre_data.SA*to_deg)
hold on
plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

%% ---FX0: fitting in pure conditions (gamma = 0, Fz = 220N)
% choose the range with: side slip = 0, camber angle = 0, vertical
% load = Fz = 220N (obv within the higher pressure dataset)

[TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );

figure('Name','FX0: pure conditions range', 'NumberTitle', 3)
plot_selected_data(TData0);

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}

FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TData0.SL));
ones_vec  = ones(size(TData0.SL));

% Initial guess (import from initialization of tyre_coeffs)
[FX0_guess,~] = MF96_FX0_vec(TData0.SL,zeros_vec , zeros_vec, tyre_coeffs_pl.FZ0*ones_vec, tyre_coeffs_pl);

% Check pure longitudinal force guess
figure('Name','FX0: guess', 'NumberTitle', 4)
plot(TData0.SL,TData0.FX,'.')
hold on
plot(TData0.SL,FX0_guess,'.')
hold off
legend({'Raw data','Guess FX0'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

% Guess values for parameters to be optimised
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
P0 = [  1,   2,   1,  0,   0,   1,   0]; 

% Limits for parameters to be optimised
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
lb = [1,   0.1,   0,   0,  -10,    0,   -10];
ub = [2,    4,   1,   1,   10,   100,  10];


KAPPA_vec = TData0.SL;
FX_vec    = TData0.FX;

% Vector for plotting: Side slip vector from -12.5° to 12.5°
SL_vec = -0.3:0.001:0.3;

% Minimization of the residual
[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fx(P,FX_vec, KAPPA_vec,0,FZ0, tyre_coeffs_pl),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs_pl.pCx1 = P_fz_nom(1) ;
tyre_coeffs_pl.pDx1 = P_fz_nom(2) ;  
tyre_coeffs_pl.pEx1 = P_fz_nom(3) ;
tyre_coeffs_pl.pEx4 = P_fz_nom(4) ;
tyre_coeffs_pl.pHx1 = P_fz_nom(5) ; 
tyre_coeffs_pl.pKx1 = P_fz_nom(6) ;
tyre_coeffs_pl.pVx1 = P_fz_nom(7) ;

res_FX0 = resid_pure_Fx(P_fz_nom , FX_vec(1), SL_vec(1), 0 , mean(TData0.FZ), tyre_coeffs_pl);


% Plot of the optimized solution
[FX0_fz_nom_vec,~] = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs_pl);

% Result of the fitting FY0 in the pure conditions
figure('Name','FX0: fitted in pure conditions','NumberTitle', 5)
plot(TData0.SL,TData0.FX,'o')
hold on
plot(SL_vec,FX0_fz_nom_vec,'-','LineWidth',2)
legend({'Raw data','Fitted FY0'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

%% ---FX0(Fz): fitting with variable Fz
% extract data with variable load, side slip equal to 0, camber angle equal to 0
[TDataDFz, ~] = intersect_table_data( SA_0, GAMMA_0 );

% extract subsets with the same vertical load
idx.FZ_220_dFz  = 220-FZ_tol < TDataDFz.FZ & TDataDFz.FZ < 220+FZ_tol;
idx.FZ_440_dFz  = 440-FZ_tol < TDataDFz.FZ & TDataDFz.FZ < 440+FZ_tol;
idx.FZ_700_dFz  = 700-FZ_tol < TDataDFz.FZ & TDataDFz.FZ < 700+FZ_tol;
idx.FZ_900_dFz  = 900-FZ_tol < TDataDFz.FZ & TDataDFz.FZ < 900+FZ_tol;
idx.FZ_1120_dFz = 1120-FZ_tol < TDataDFz.FZ & TDataDFz.FZ < 1120+FZ_tol;
FZ_220_dFz  = TDataDFz( idx.FZ_220_dFz, : );
FZ_700_dFz  = TDataDFz( idx.FZ_700_dFz, : );
FZ_900_dFz  = TDataDFz( idx.FZ_900_dFz, : );
FZ_1120_dFz = TDataDFz( idx.FZ_1120_dFz, : );

zeros_vec = zeros(size(TDataDFz.SL));
ones_vec  = ones(size(TDataDFz.SL));

% Guess values for parameters to be optimised
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0_dFz = [  0,   0,   0,  0,   0,   0,   0]; 

% Limits for parameters to be optimised
lb = [];
ub = [];

KAPPA_vec = TDataDFz.SL;
FX_vec    = TDataDFz.FX;
FZ_vec    = TDataDFz.FZ;

% check guess
[FX0_dfz_vec,~] = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                           tyre_coeffs_pl.FZ0*ones(size(SL_vec)),tyre_coeffs_pl);

figure('Name','FX0(Fz): guess', 'NumberTitle', 6)
plot(KAPPA_vec,FX_vec,'.')
hold on
plot(SL_vec,FX0_dfz_vec,'.')
legend({'Raw data variable load','Guess FX0'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}(Fz)$ [N]')

% Residual minimization
[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varFz(P,FX_vec, KAPPA_vec,0,FZ_vec, tyre_coeffs_pl),...
                               P0_dFz,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs_pl.pDx2 = P_dfz(1) ;
tyre_coeffs_pl.pEx2 = P_dfz(2) ;  
tyre_coeffs_pl.pEx3 = P_dfz(3) ;
tyre_coeffs_pl.pHx2 = P_dfz(4) ;
tyre_coeffs_pl.pKx2 = P_dfz(5) ; 
tyre_coeffs_pl.pKx3 = P_dfz(6) ;
tyre_coeffs_pl.pVx2 = P_dfz(7) ; 

res_FX0_dfz_vec = resid_pure_Fx_varFz(P_dfz,FX_vec,SL_vec,0 , FZ_vec,tyre_coeffs_pl);

tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));

[FX0_fz_var_vec1, Kxk_fz_var_vec1] = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs_pl);
[FX0_fz_var_vec2, Kxk_fz_var_vec2] = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs_pl);
[FX0_fz_var_vec3, Kxk_fz_var_vec3] = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs_pl);
[FX0_fz_var_vec4, Kxk_fz_var_vec4] = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs_pl);


figure('Name','FX0(Fz): fitted with variable Fz','NumberTitle', 7)
hold on
plot(FZ_220_dFz.SL,FZ_220_dFz.FX,'.','Color',[0.8500 0.3250 0.0980])
plot(FZ_700_dFz.SL,FZ_700_dFz.FX,'.','Color',[0.6350 0.0780 0.1840])
plot(FZ_900_dFz.SL,FZ_900_dFz.FX,'.','Color',[0.4940 0.1840 0.5560])
plot(FZ_1120_dFz.SL,FZ_1120_dFz.FX,'.','Color',[0.3010 0.7450 0.9330])
plot(SL_vec,FX0_fz_var_vec1,'-','LineWidth',2,'Color',[1 0.7 0.34])
plot(SL_vec,FX0_fz_var_vec2,'-','LineWidth',2,'Color','r')
plot(SL_vec,FX0_fz_var_vec3,'-','LineWidth',2,'Color','m')
plot(SL_vec,FX0_fz_var_vec4,'-','LineWidth',2,'Color','c')
hold off
legend({'Raw with $Fz=220N$','Raw with $Fz=700N$','Raw with $Fz=900N$','Raw with $Fz=1120N$', '$Fx(Fz_{220})$','$Fx(Fz_{700})$','$Fx(Fz_{900})$','$Fx(Fz_{1120})$'}, 'Location','eastoutside');
xlabel('$\kappa$ [.]')
ylabel('$F_{x0}(Fz)$ [N]')

figure('Name','Kxk(Fz): cornering stiffness as function of Fz','NumberTitle', 8)
hold on
plot(mean(FZ_220_dFz.FZ),Kxk_fz_var_vec1(1),'+')
plot(mean(FZ_700_dFz.FZ),Kxk_fz_var_vec2(1),'+')
plot(mean(FZ_900_dFz.FZ),Kxk_fz_var_vec3(1),'+')
plot(mean(FZ_1120_dFz.FZ),Kxk_fz_var_vec4(1),'+')
hold off
legend({'Kxk($Fz_{220}$)','Kxk($Fz_{700}$)','Kxk($Fz_{900}$)','Kxk($Fz_{1120}$)'}, 'Location','eastoutside');
xlabel('$Fz$ [N]')
ylabel('$K_{kx}(Fz)$ [-]')

Calfa_vec1_x = MF96_CorneringStiffness_x(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220_dFz.FZ)*tmp_ones,tyre_coeffs_pl);
Calfa_vec2_x = MF96_CorneringStiffness_x(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700_dFz.FZ)*tmp_ones,tyre_coeffs_pl);
Calfa_vec3_x = MF96_CorneringStiffness_x(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900_dFz.FZ)*tmp_ones,tyre_coeffs_pl);
Calfa_vec4_x = MF96_CorneringStiffness_x(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120_dFz.FZ)*tmp_ones,tyre_coeffs_pl);

figure('Name','Kxk(kappa): cornering stiffness as function of kappa','NumberTitle', 9)
hold on
plot(SL_vec,Calfa_vec1_x,'-','LineWidth',2)
plot(SL_vec,Calfa_vec2_x,'-','LineWidth',2)
plot(SL_vec,Calfa_vec3_x,'-','LineWidth',2)
plot(SL_vec,Calfa_vec4_x,'-','LineWidth',2)
hold off
legend({'Kxk($Fz_{220}$)','Kxk($Fz_{700}$)','Kxk($Fz_{900}$)','Kxk($Fz_{1120}$)'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$K_{kx}(Fz)$ [-]')

%% ---FX0(gamma): fitting with variable camber(gamma)
% extract data with the same vertical load (Fz = 220N)
[TDataGamma, ~] = intersect_table_data( SA_0, FZ_220 );

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
P0_dgamma = [0]; 

% Limits for parameters to be optimised
lb = [];
ub = [];

zeros_vec = zeros(size(TDataGamma.SL));
ones_vec  = ones(size(TDataGamma.SL));

KAPPA_vec = TDataGamma.SL;
GAMMA_vec = TDataGamma.IA; 
FX_vec    = TDataGamma.FX;
FZ_vec    = TDataGamma.FZ;

[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs_pl.FZ0, tyre_coeffs_pl),...
                               P0_dgamma,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs_pl.pDx3 = P_varGamma(1) ; % 1

[FX0_varGamma_vec,~] = MF96_FX0_vec(KAPPA_vec,zeros_vec , GAMMA_vec, tyre_coeffs_pl.FZ0*ones_vec,tyre_coeffs_pl);

figure('Name','FX0(Gamma): fitted with variable camber','NumberTitle', 10)
plot(KAPPA_vec,TDataGamma.FX,'o')
hold on
plot(KAPPA_vec,FX0_varGamma_vec,'-')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

% Calculate the residuals with the optimal solution found above
res_Fx0_varGamma  = resid_pure_Fx_varGamma(P_varGamma,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs_pl.FZ0, tyre_coeffs_pl);

%% ---------------last figure FX0---------------
last_fig_FX0 = 10;

%% --Pure lateral force FY0: dataset import

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

%% ---Dataset for pure lateral: plot

figure ('Name','FY0: entire raw dataset', 'NumberTitle', 1 + last_fig_FX0)
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

%% ---Higher pressure dataset for pure lateral force: table selection and plot
% consider the high pressure region of original dataset (more stable one)

vec_samples_pl = 1:1:length(smpl_range_pl);

tyre_data_pl = table();
% store raw data in table
tyre_data_pl.SL =  SL(smpl_range_pl);
tyre_data_pl.SA = -SA(smpl_range_pl)*to_rad;    % SAE -> Adapted SAE
tyre_data_pl.FZ = -FZ(smpl_range_pl);           % SAE -> Adapted SAE
tyre_data_pl.FX =  FX(smpl_range_pl);
tyre_data_pl.FY =  FY(smpl_range_pl);   
tyre_data_pl.MZ =  MZ(smpl_range_pl);
tyre_data_pl.IA =  IA(smpl_range_pl)*to_rad;

% Extract points at constant camber angle

% Test data done at: 
%  - 0 deg
%  - 1 deg
%  - 2 deg
%  - 3 deg
%  - 4 deg
% in the following order: (0 2 4 1 3)*2

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
%  - 100lbf (100*0.453592*9.81 =  445N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )
% in the following order: (200 150 50 250 100)*2

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
figure('Name','FY0: higher pressure dataset with regions', 'NumberTitle', 2 + last_fig_FX0)
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

%% ---FY0: fitting in pure conditions (gamma = 0, Fz = 220N)
% choose the range with: longitudinal slip = 0, camber angle = 0, vertical
% load = Fz = 220N (obv within the higher pressure dataset)

[TData0_pl, ~] = intersect_table_data( GAMMA_0_pl, FZ_220_pl );

figure('Name','FY0: pure conditions range', 'NumberTitle', 3 + last_fig_FX0)
plot_selected_data(TData0_pl);

% Fit the coefficients {pCy1, pDy1, pEy1, pHy1, pKy1, pKy2, pVy1}

zeros_vec_pl = zeros(size(TData0_pl.SA));
ones_vec_pl  = ones(size(TData0_pl.SA));

% Initial guess (import from initialization of tyre_coeffs)
[FY0_guess,~] = MF96_FY0_vec(zeros_vec_pl, TData0_pl.SA , zeros_vec_pl, tyre_coeffs_pl.FZ0*ones_vec_pl, tyre_coeffs_pl);

% Check lateral pure force guess
figure('Name','FY0: guess', 'NumberTitle', 4 + last_fig_FX0)
plot(TData0_pl.SA*to_deg,TData0_pl.FY,'.')
hold on
plot(TData0_pl.SA*to_deg,FY0_guess,'.')
hold off
legend({'Raw data','Guess FY0'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}$ [N]')

% Guess values for parameters to be optimised
%       [pCy1 pDy1 pEy1 pHy1  pKy1  pKy2  pVy1]
P0_pl = [ 1.5, 2.4, 0.34, -0.48e-2, -0.29e2, -0.1e1, 0.47e-1 ]; 

% Limits for parameters to be optimised
lb_pl = [ ];
ub_pl = [ ];

ALPHA_vec = TData0_pl.SA;
FY_vec    = TData0_pl.FY;

% Vector for plotting: Side slip vector from -12.5° to 12.5°
SA_vec = (-12.5*to_rad):0.001:(12.5*to_rad);

% Minimization of the residual
[P_fz_nom_pl,fval,exitflag] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,mean(TData0_pl.FZ), tyre_coeffs_pl),...
                               P0_pl,[],[],[],[],lb_pl,ub_pl);

% Update tyre data with new optimal values                            
tyre_coeffs_pl.pCy1 = P_fz_nom_pl(1) ;
tyre_coeffs_pl.pDy1 = P_fz_nom_pl(2) ;  
tyre_coeffs_pl.pEy1 = P_fz_nom_pl(3) ;
tyre_coeffs_pl.pHy1 = P_fz_nom_pl(4) ;
tyre_coeffs_pl.pKy1 = P_fz_nom_pl(5) ; 
tyre_coeffs_pl.pKy2 = P_fz_nom_pl(6) ;
tyre_coeffs_pl.pVy1 = P_fz_nom_pl(7) ;

res_FY0 = resid_pure_Fy(P_fz_nom_pl , FY_vec, SA_vec, 0 , mean(TData0_pl.FZ), tyre_coeffs_pl);


% Plot of the optimized solution
[FY0_fz_nom_vec,~] = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec , zeros(size(SA_vec)), ...
                              mean(TData0_pl.FZ).*ones(size(SA_vec)),tyre_coeffs_pl);

% Result of the fitting FY0 in the pure conditions
figure('Name','FY0: fitted in pure conditions','NumberTitle', 5 + last_fig_FX0)
plot(TData0_pl.SA*to_deg,TData0_pl.FY,'o')
hold on
plot(SA_vec*to_deg,FY0_fz_nom_vec,'-','LineWidth',2)
legend({'Raw data','Fitted FY0'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}$ [N]')

%% ---FY0(Fz): fitting with variable Fz
% extract data with variable load and camber angle equal to 0
TDataDFz_pl = GAMMA_0_pl;

% figure('Name','Variable Fz range', 'NumberTitle', 6 + last_fig_FX0)
% plot_selected_data(TDataDFz);

smpl_range_pl_dFz = size(TDataDFz_pl);
vec_samples_pl_dFz = 1:1:smpl_range_pl_dFz;

% Extract points at constant vertical load
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
figure('Name','FY0(Fz): considered dataset', 'NumberTitle', 6 + last_fig_FX0)
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

% Guess values for parameters to be optimised
%    [pDy2 pEy2 pHy2 pVy2] 
P0_pl_dFz =[ 0, 0, 0, 0 ];   

% Limits for parameters to be optimised
%    [pDy2 pEy2 pHy2 pVy2] 
lb_dFz = [ ];
ub_dFz = [ ];

ALPHA_vec_dFz = TDataDFz_pl.SA;
FY_vec_dFz    = TDataDFz_pl.FY;
FZ_vec_dFz    = TDataDFz_pl.FZ;

% check guess for variable load
[FY0_dfz_vec,~] = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              mean(FZ_220_pl_dFz.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);

figure('Name','FY0(Fz): guess', 'NumberTitle', 7 + last_fig_FX0)
plot(ALPHA_vec_dFz*to_deg,FY_vec_dFz,'.')
hold on
plot(SA_vec*to_deg,FY0_dfz_vec,'.')
legend({'Raw data variable load','Guess FY0'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}(Fz)$ [N]')

% Resitual minimization
[P_dfz_pl,fval,exitflag] = fmincon(@(P_pl)resid_pure_Fy_varFz(P_pl,FY_vec_dFz, ALPHA_vec_dFz,0,FZ_vec_dFz, tyre_coeffs_pl),...
                               P0_pl_dFz,[],[],[],[],lb_dFz,ub_dFz);

% Change tyre data with new optimal values                             
tyre_coeffs_pl.pDy2 = P_dfz_pl(1);
tyre_coeffs_pl.pEy2 = P_dfz_pl(2);
tyre_coeffs_pl.pHy2 = P_dfz_pl(3);
tyre_coeffs_pl.pVy2 = P_dfz_pl(4);

res_FY0_dfz = resid_pure_Fy_varFz(P_dfz_pl , FY_vec_dFz,SA_vec, 0 , FZ_vec_dFz, tyre_coeffs_pl);

tmp_zeros_dFz = zeros(size(SA_vec));
tmp_ones_dFz = ones(size(SA_vec));

[FY0_fz_var_vec1, Kya_fz_var_vec1] = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_220_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
[FY0_fz_var_vec2, Kya_fz_var_vec2] = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_440_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
[FY0_fz_var_vec3, Kya_fz_var_vec3] = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_700_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
[FY0_fz_var_vec4, Kya_fz_var_vec4] = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_900_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
[FY0_fz_var_vec5, Kya_fz_var_vec5] = MF96_FY0_vec(tmp_zeros_dFz, SA_vec ,tmp_zeros_dFz, mean(FZ_1120_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);


figure('Name','FY0(Fz): fitted with variable Fz','NumberTitle', 8 + last_fig_FX0)
hold on
plot(FZ_220_pl_dFz.SA*to_deg,FZ_220_pl_dFz.FY,'o', 'Color', [1 0.4 0.4])
plot(FZ_440_pl_dFz.SA*to_deg,FZ_440_pl_dFz.FY,'o', 'Color', [0.6 0.8 1])
plot(FZ_700_pl_dFz.SA*to_deg,FZ_700_pl_dFz.FY,'o', 'Color', [0.5 1 0.5])
plot(FZ_900_pl_dFz.SA*to_deg,FZ_900_pl_dFz.FY,'o', 'Color', [0.8 0.6 1])
plot(FZ_1120_pl_dFz.SA*to_deg,FZ_1120_pl_dFz.FY,'o', 'Color', [0.6 1 1])
plot(SA_vec*to_deg,FY0_fz_var_vec1,'.','LineWidth',2, 'Color', 'r')
plot(SA_vec*to_deg,FY0_fz_var_vec2,'.','LineWidth',2, 'Color', 'b')
plot(SA_vec*to_deg,FY0_fz_var_vec3,'.','LineWidth',2, 'Color', 'g')
plot(SA_vec*to_deg,FY0_fz_var_vec4,'.','LineWidth',2, 'Color', 'm')
plot(SA_vec*to_deg,FY0_fz_var_vec5,'.','LineWidth',2, 'Color', 'c')
hold off
legend({'Raw with $Fz=220N$','Raw with $Fz=440N$','Raw with $Fz=700N$','Raw with $Fz=900N$','Raw with $Fz=1120N$', '$Fy(Fz_{220})$','$Fy(Fz_{440})$','$Fy(Fz_{700})$','$Fy(Fz_{900})$','$Fy(Fz_{1120})$'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}(Fz)$ [N]')

% Stiffness
figure('Name','Kya(Fz): cornering stiffness as function of Fz','NumberTitle', 9 + last_fig_FX0)
hold on
plot(mean(FZ_220_pl_dFz.FZ),Kya_fz_var_vec1(1),'+','MarkerSize',10)
plot(mean(FZ_440_pl_dFz.FZ),Kya_fz_var_vec2(1),'+','MarkerSize',10)
plot(mean(FZ_700_pl_dFz.FZ),Kya_fz_var_vec3(1),'+','MarkerSize',10)
plot(mean(FZ_900_pl_dFz.FZ),Kya_fz_var_vec4(1),'+','MarkerSize',10)
plot(mean(FZ_1120_pl_dFz.FZ),Kya_fz_var_vec5(1),'+','MarkerSize',10)
hold off
legend({'Kya($Fz_{220}$)','Kya($Fz_{440}$)','Kya($Fz_{700}$)','Kya($Fz_{900}$)','Kya($Fz_{1120}$)'}, 'Location','eastoutside');
xlabel('$Fz$ [N]')
ylabel('$K_{ya}(Fz)$')

Calfa_vec1_y = MF96_CorneringStiffness_y(tmp_zeros_dFz,SA_vec ,tmp_zeros_dFz, mean(FZ_220_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
Calfa_vec2_y = MF96_CorneringStiffness_y(tmp_zeros_dFz,SA_vec ,tmp_zeros_dFz, mean(FZ_440_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
Calfa_vec3_y = MF96_CorneringStiffness_y(tmp_zeros_dFz,SA_vec ,tmp_zeros_dFz, mean(FZ_700_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
Calfa_vec4_y = MF96_CorneringStiffness_y(tmp_zeros_dFz,SA_vec ,tmp_zeros_dFz, mean(FZ_900_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);
Calfa_vec5_y = MF96_CorneringStiffness_y(tmp_zeros_dFz,SA_vec ,tmp_zeros_dFz, mean(FZ_1120_pl_dFz.FZ)*tmp_ones_dFz,tyre_coeffs_pl);

figure('Name','Kya(alpha): cornering stiffness as function of alpha','NumberTitle', 10 + last_fig_FX0)
hold on
plot(SA_vec*to_deg,Calfa_vec1_y,'-','LineWidth',2)
plot(SA_vec*to_deg,Calfa_vec2_y,'-','LineWidth',2)
plot(SA_vec*to_deg,Calfa_vec3_y,'-','LineWidth',2)
plot(SA_vec*to_deg,Calfa_vec4_y,'-','LineWidth',2)
plot(SA_vec*to_deg,Calfa_vec5_y,'-','LineWidth',2)
hold off
legend({'Kya($Fz_{220}$)','Kya($Fz_{440}$)','Kya($Fz_{700}$)','Kya($Fz_{900}$)','Kya($Fz_{1120}$)'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$K_{kx}(Fz)$ [-]')

%% ---FY0(gamma): fitting with variable camber(gamma)
% extract data with the same vertical load (Fz = 220N) 
TDataGamma_pl = FZ_220_pl;

smpl_range_pl_dgamma = size(TDataGamma_pl);
vec_samples_pl_dgamma = 1:1:smpl_range_pl_dgamma;

% Extract points at constant camber
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
figure('Name','FY0(gamma): considered dataset', 'NumberTitle', 11 + last_fig_FX0)
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
%   [pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4]
P0_pl_dgamma = [ 0 , 0 , 0 , 0 , 0 , 0 , 0 ];

% Limits for parameters to be optimised
lb_dgamma = [ ];
ub_dgamma = [ ];

zeros_vec_dgamma = zeros(size(TDataGamma_pl.IA));
ones_vec_dgamma  = ones(size(TDataGamma_pl.IA));

ALPHA_vec_dgamma = TDataGamma_pl.SA; 
GAMMA_vec_dgamma = TDataGamma_pl.IA; 
FY_vec_dgamma    = TDataGamma_pl.FY;
FZ_vec_dgamma    = TDataGamma_pl.FZ;

[FY0_varGamma_vec,~] = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec , GAMMA_vec_dgamma, tyre_coeffs_pl.FZ0*ones(size(SA_vec)),tyre_coeffs_pl);

figure('Name','FY0(gamma): guess', 'NumberTitle', 12 + last_fig_FX0)
plot(ALPHA_vec_dgamma,TDataGamma_pl.FY,'o')
hold on
plot(SA_vec,FY0_varGamma_vec,'.','MarkerSize',5)
legend({'Raw data variable camber','Guess FY0(gamma)'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}(\gamma)$ [N]')

% Residual minimization
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

tmp_zeros_dgamma = zeros(size(SA_vec));
tmp_ones_dgamma = ones(size(SA_vec));

[FY0_gamma_var_vec1,~] = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_0_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
[FY0_gamma_var_vec2,~] = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_1_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
[FY0_gamma_var_vec3,~] = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_2_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
[FY0_gamma_var_vec4,~] = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_3_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);
[FY0_gamma_var_vec5,~] = MF96_FY0_vec(tmp_zeros_dgamma, SA_vec ,mean(GAMMA_4_dgamma.IA)*tmp_ones_dgamma, mean(TDataGamma_pl.FZ)*tmp_ones_dgamma,tyre_coeffs_pl);


figure('Name','FY0(gamma): fitted with variable camber','NumberTitle', 13 + last_fig_FX0)
hold on
plot(GAMMA_0_dgamma.SA*to_deg,GAMMA_0_dgamma.FY,'-','MarkerSize',5, 'Color', '#0072BD') %'MarkerEdgeColor','y',
plot(GAMMA_1_dgamma.SA*to_deg,GAMMA_1_dgamma.FY,'-','MarkerSize',5, 'Color', '#D95319') %'MarkerEdgeColor','c',
plot(GAMMA_2_dgamma.SA*to_deg,GAMMA_2_dgamma.FY,'-','MarkerSize',5, 'Color', '#EDB120') %'MarkerEdgeColor','m',
plot(GAMMA_3_dgamma.SA*to_deg,GAMMA_3_dgamma.FY,'-','MarkerSize',5, 'Color', '#77AC30') %'MarkerEdgeColor','b',
plot(GAMMA_4_dgamma.SA*to_deg,GAMMA_4_dgamma.FY,'-','MarkerSize',5, 'Color', '#4DBEEE') %'MarkerEdgeColor','r',
plot(SA_vec*to_deg,FY0_gamma_var_vec1,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#0072BD')
plot(SA_vec*to_deg,FY0_gamma_var_vec2,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#D95319')
plot(SA_vec*to_deg,FY0_gamma_var_vec3,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#EDB120')
plot(SA_vec*to_deg,FY0_gamma_var_vec4,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#77AC30')
plot(SA_vec*to_deg,FY0_gamma_var_vec5,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#4DBEEE')
legend({'Raw data with $\gamma = 0 deg $','Raw data with $\gamma = 1 deg $','Raw data with $\gamma = 2 deg $','Raw data with $\gamma = 3 deg $','Raw data with $\gamma = 4 deg $', 'Fy($\gamma = 0 deg$)','Fy($\gamma = 1 deg$)','Fy($\gamma = 2 deg$)','Fy($\gamma = 3 deg$)','Fy($\gamma = 4 deg$)'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}(\gamma)$ [N]')

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

%% ---------------last figure FY0---------------
last_fig_FY0 = 13 + last_fig_FX0;

%% --Pure self aligning moment MZO: same dataset (lateral)
cut_start_mz = 27760;
cut_end_mz = 54500;

% select dataset portion
smpl_range_mz = cut_start_mz:cut_end_mz;

%% ---Dataset for pure self-aligning moment MZ0: plot

figure ('Name','MZ0: entire raw dataset', 'NumberTitle', 1 + last_fig_FY0)
tiledlayout(6,1)

ax_list_mz(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start_mz cut_start_mz],y_range,'--r')
plot([cut_end_mz cut_end_mz],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list_mz(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start_mz cut_start_mz],y_range,'--r')
plot([cut_end_mz cut_end_mz],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_mz(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start_mz cut_start_mz],y_range,'--r')
plot([cut_end_mz cut_end_mz],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_mz(4) = nexttile; y_range = [min(min(MZ),0) round(max(MZ)*1.1)];
plot(MZ)
hold on
plot([cut_start_mz cut_start_mz],y_range,'--r')
plot([cut_end_mz cut_end_mz],y_range,'--r')
title('Aligning moment')
xlabel('Samples [-]')
ylabel('[Nm]')

ax_list_mz(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start_mz cut_start_mz],y_range,'--r')
plot([cut_end_mz cut_end_mz],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[Pa]')

ax_list_mz(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start_mz cut_start_mz],y_range,'--r')
plot([cut_end_mz cut_end_mz],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list_mz,'x')

%% ---Higher pressure dataset for pure self-aligning moment: table selection and plot
% consider the high pressure region of original dataset (more stable one)

vec_samples_mz = 1:1:length(smpl_range_mz);
tyre_data_mz = table(); % create empty table
% store raw data in table
tyre_data_mz.SL =  SL(smpl_range_mz);
tyre_data_mz.SA = -SA(smpl_range_mz)*to_rad;
tyre_data_mz.FZ = -FZ(smpl_range_mz); 
tyre_data_mz.FX =  FX(smpl_range_mz);
tyre_data_mz.FY =  FY(smpl_range_mz);
tyre_data_mz.MZ =  MZ(smpl_range_mz);
tyre_data_mz.IA =  IA(smpl_range_mz)*to_rad;

% Extract points at constant camber angle

% Test data done at: 
%  - 0 deg
%  - 1 deg
%  - 2 deg
%  - 3 deg
%  - 4 deg
% in the following order: (0 2 4 1 3)*2

GAMMA_tol_mz = 0.05*to_rad;
idx_mz.GAMMA_0 = 0.0*to_rad-GAMMA_tol_mz < tyre_data_mz.IA & tyre_data_mz.IA < 0.0*to_rad+GAMMA_tol_mz;
idx_mz.GAMMA_1 = 1.0*to_rad-GAMMA_tol_mz < tyre_data_mz.IA & tyre_data_mz.IA < 1.0*to_rad+GAMMA_tol_mz;
idx_mz.GAMMA_2 = 2.0*to_rad-GAMMA_tol_mz < tyre_data_mz.IA & tyre_data_mz.IA < 2.0*to_rad+GAMMA_tol_mz;
idx_mz.GAMMA_3 = 3.0*to_rad-GAMMA_tol_mz < tyre_data_mz.IA & tyre_data_mz.IA < 3.0*to_rad+GAMMA_tol_mz;
idx_mz.GAMMA_4 = 4.0*to_rad-GAMMA_tol_mz < tyre_data_mz.IA & tyre_data_mz.IA < 4.0*to_rad+GAMMA_tol_mz;
GAMMA_0_mz  = tyre_data_mz( idx_mz.GAMMA_0, : );
GAMMA_1_mz  = tyre_data_mz( idx_mz.GAMMA_1, : );
GAMMA_2_mz  = tyre_data_mz( idx_mz.GAMMA_2, : );
GAMMA_3_mz  = tyre_data_mz( idx_mz.GAMMA_3, : );
GAMMA_4_mz  = tyre_data_mz( idx_mz.GAMMA_4, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol_mz = 100;
idx_mz.FZ_220  = 220-FZ_tol_mz < tyre_data_mz.FZ & tyre_data_mz.FZ < 220+FZ_tol_mz;
idx_mz.FZ_440  = 440-FZ_tol_mz < tyre_data_mz.FZ & tyre_data_mz.FZ < 440+FZ_tol_mz;
idx_mz.FZ_700  = 700-FZ_tol_mz < tyre_data_mz.FZ & tyre_data_mz.FZ < 700+FZ_tol_mz;
idx_mz.FZ_900  = 900-FZ_tol_mz < tyre_data_mz.FZ & tyre_data_mz.FZ < 900+FZ_tol_mz;
idx_mz.FZ_1120 = 1120-FZ_tol_mz < tyre_data_mz.FZ & tyre_data_mz.FZ < 1120+FZ_tol_mz;
FZ_220_mz  = tyre_data_mz( idx_mz.FZ_220, : );
FZ_440_mz  = tyre_data_mz( idx_mz.FZ_440, : );
FZ_700_mz  = tyre_data_mz( idx_mz.FZ_700, : );
FZ_900_mz  = tyre_data_mz( idx_mz.FZ_900, : );
FZ_1120_mz = tyre_data_mz( idx_mz.FZ_1120, : );

% Plot
figure('Name','MZ0: higher pressure dataset with regions', 'NumberTitle', 2 + last_fig_FY0)
tiledlayout(2,1)

ax_list_mz_2(1) = nexttile;
plot(tyre_data_mz.IA*to_deg)
hold on
plot(vec_samples_mz(idx_mz.GAMMA_0),GAMMA_0_mz.IA*to_deg,'.');
plot(vec_samples_mz(idx_mz.GAMMA_1),GAMMA_1_mz.IA*to_deg,'.');
plot(vec_samples_mz(idx_mz.GAMMA_2),GAMMA_2_mz.IA*to_deg,'.');
plot(vec_samples_mz(idx_mz.GAMMA_3),GAMMA_3_mz.IA*to_deg,'.');
plot(vec_samples_mz(idx_mz.GAMMA_4),GAMMA_4_mz.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_mz_2(2) = nexttile;
plot(tyre_data_mz.FZ)
hold on
plot(vec_samples_mz(idx_mz.FZ_220),FZ_220_mz.FZ,'.');
plot(vec_samples_mz(idx_mz.FZ_440),FZ_440_mz.FZ,'.');
plot(vec_samples_mz(idx_mz.FZ_700),FZ_700_mz.FZ,'.');
plot(vec_samples_mz(idx_mz.FZ_900),FZ_900_mz.FZ,'.');
plot(vec_samples_mz(idx_mz.FZ_1120),FZ_1120_mz.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

linkaxes(ax_list_mz_2,'x')

%% ---MZ0: fitting in pure conditions (gamma = 0, Fz = 220N)
% choose the range with: longitudinal slip = 0, camber angle = 0, vertical
% load = Fz = 220N (obv within the higher pressure dataset)

[TData0_mz, ~] = intersect_table_data( GAMMA_0_mz, FZ_220_mz );

figure('Name','MZ0: pure conditions range', 'NumberTitle', 3 + last_fig_FX0)
plot_selected_data(TData0_mz);

% Fit the coefficients {qBz1, qBz9, qBz10, qCz1, qDz1, qDz6, qEz1, qEz4, qHz1}

zeros_vec_mz = zeros(size(TData0_mz.MZ));
ones_vec_mz  = ones(size(TData0_mz.MZ));

% Initial guess (import from initialization of tyre_coeffs)
MZ0_guess = MF96_MZ0_vec(zeros_vec_mz, TData0_mz.SA, zeros_vec_mz, tyre_coeffs_pl.FZ0*ones_vec_mz, tyre_coeffs_pl);


% Check lateral pure force guess
figure('Name','MZ0: guess', 'NumberTitle', 4 + last_fig_FY0)
plot(TData0_mz.SA*to_deg,TData0_mz.MZ,'.')
hold on
plot(TData0_mz.SA*to_deg,MZ0_guess,'.')
hold off
legend({'Raw data','Guess MZ0'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}$ [Nm]')

% Guess values for parameters to be optimised
%         {qBz1, qBz9, qBz10, qCz1,  qDz1,  qDz6, qEz1,  qEz4, qHz1}
P0_mz = [  3.11, 1.87, 12.10, 2.82, -0.96, -0.01, 2.90, -0.33, 0.48];

% Limits for parameters to be optimised
lb_mz_dgamma = [ ];
ub_mz_dgamma = [ ];


ALPHA_vec_mz = TData0_mz.SA;
MZ_vec_mz    = TData0_mz.MZ;


% Residual minimization
[P_fz_nom_mz,fval,exitflag] = fmincon(@(P)resid_pure_Mz(P, MZ_vec_mz, ALPHA_vec_mz, 0, mean(TData0_mz.FZ), tyre_coeffs_pl),...
                               P0_mz,[],[],[],[],lb_mz_dgamma,ub_mz_dgamma);

% Update tyre data with new optimal values                             
tyre_coeffs_pl.qBz1 = P_fz_nom_mz(1); 
tyre_coeffs_pl.qBz9 = P_fz_nom_mz(2);
tyre_coeffs_pl.qBz10 = P_fz_nom_mz(3);
tyre_coeffs_pl.qCz1 = P_fz_nom_mz(4);
tyre_coeffs_pl.qDz1 = P_fz_nom_mz(5);
tyre_coeffs_pl.qDz6 = P_fz_nom_mz(6);
tyre_coeffs_pl.qEz1 = P_fz_nom_mz(7);
tyre_coeffs_pl.qEz4 = P_fz_nom_mz(8);
tyre_coeffs_pl.qHz1 = P_fz_nom_mz(9);

% Plot of the optimized solution
MZ0_fz_nom_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              mean(TData0_mz.FZ).*ones(size(SA_vec)),tyre_coeffs_pl);


% Result of the fitting MZ0 in the pure conditions
figure('Name','MZ0: fitted in pure conditions','NumberTitle', 5 + last_fig_FY0)
plot(TData0_mz.SA*to_deg,TData0_mz.MZ,'o')
hold on
plot(SA_vec*to_deg,MZ0_fz_nom_vec,'-','LineWidth',2)
legend({'Raw data','Fitted FY0'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}$ [Nm]')

%% ---MZ0(Fz): fitting with variable Fz
% extract data with variable load and camber angle equal to 0

TDataDFz_mz = GAMMA_0_mz;

smpl_range_mz_dFz = size(TDataDFz_mz);

% Extract points at constant vertical load
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
FZ_1120_mz_dFz = TDataDFz_mz( idx_mz_dFz.FZ_1120, : );


% Fit the coeffs {qBz2, qBz3, qDz2, qDz7, qEz2, qEz3, qHz2}

zeros_vec_mz = zeros(size(TDataDFz_mz.SA));
ones_vec_mz  = ones(size(TDataDFz_mz.SA));


% Guess values for parameters to be optimised
% [qBz2, qBz3, qDz2, qDz7, qEz2, qEz3, qHz2]
P0_mz_dFz = [ 0,  0,  0,  0,  0,  0,   0 ]; 

% Limits for parameters to be optimised
lb_mz_dFz = [];
ub_mz_dFz = [];


ALPHA_vec_mz_dFz = TDataDFz_mz.SA;
Mz_vec_dFz    = TDataDFz_mz.MZ;
FZ_vec_dFz    = TDataDFz_mz.FZ;

% check guess
MZ0_dfz_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                           FZ_vec_dFz, tyre_coeffs_pl);

figure('Name','MZ0(Fz): guess', 'NumberTitle', 7 + last_fig_FY0)
plot(ALPHA_vec_mz_dFz*to_deg,Mz_vec_dFz,'.')
hold on
plot(SA_vec*to_deg,MZ0_dfz_vec,'.')
legend({'Raw data variable load','Guess MZ0'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}(Fz)$ [Nm]')

% Residual minimization
[P_dfz_mz,fval,exitflag] = fmincon(@(P)resid_pure_Mz_varFz(P, Mz_vec_dFz, ALPHA_vec_mz_dFz, 0, FZ_vec_dFz, tyre_coeffs_pl),...
                               P0_mz_dFz,[],[],[],[],lb_mz_dFz,ub_mz_dFz);

% Change tyre data with new optimal values                             
tyre_coeffs_pl.qBz2 = P_dfz_mz(1); 
tyre_coeffs_pl.qBz3 = P_dfz_mz(2);
tyre_coeffs_pl.qDz2 = P_dfz_mz(3);
tyre_coeffs_pl.qDz7 = P_dfz_mz(4);
tyre_coeffs_pl.qEz2 = P_dfz_mz(5);
tyre_coeffs_pl.qEz3 = P_dfz_mz(6);
tyre_coeffs_pl.qHz2 = P_dfz_mz(7);


res_MZ0_dfz_vec = resid_pure_Mz_varFz(P_dfz_mz, Mz_vec_dFz, SA_vec, 0 , FZ_vec_dFz, tyre_coeffs_pl);

tmp_zeros_mz_dFz = zeros(size(SA_vec));
tmp_ones_mz_dFz = ones(size(SA_vec));

MZ0_fz_var_vec1 = MF96_MZ0_vec(tmp_zeros_mz_dFz, SA_vec, tmp_zeros_mz_dFz, mean(FZ_220_mz.FZ)*tmp_ones_mz_dFz,tyre_coeffs_pl);
MZ0_fz_var_vec2 = MF96_MZ0_vec(tmp_zeros_mz_dFz, SA_vec, tmp_zeros_mz_dFz, mean(FZ_440_mz.FZ)*tmp_ones_mz_dFz,tyre_coeffs_pl);
MZ0_fz_var_vec3 = MF96_MZ0_vec(tmp_zeros_mz_dFz, SA_vec, tmp_zeros_mz_dFz, mean(FZ_700_mz.FZ)*tmp_ones_mz_dFz,tyre_coeffs_pl);
MZ0_fz_var_vec4 = MF96_MZ0_vec(tmp_zeros_mz_dFz, SA_vec, tmp_zeros_mz_dFz, mean(FZ_900_mz.FZ)*tmp_ones_mz_dFz,tyre_coeffs_pl);
MZ0_fz_var_vec5 = MF96_MZ0_vec(tmp_zeros_mz_dFz, SA_vec, tmp_zeros_mz_dFz, mean(FZ_1120_mz.FZ)*tmp_ones_mz_dFz,tyre_coeffs_pl);


figure('Name','Mz0(Fz0)')
plot(TDataDFz_mz.SA*to_deg,TDataDFz_mz.MZ,'o')
hold on

plot(SA_vec*to_deg,MZ0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec*to_deg,MZ0_fz_var_vec5,'-','LineWidth',2)

xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}$ [N]')

figure('Name','MZ0(Fz): fitted with variable Fz','NumberTitle', 8 + last_fig_FY0)
hold on
plot(FZ_220_mz_dFz.SA*to_deg,FZ_220_mz_dFz.MZ,'o', 'Color', [1 0.4 0.4])
plot(FZ_440_mz_dFz.SA*to_deg,FZ_440_mz_dFz.MZ,'o', 'Color', [0.6 0.8 1])
plot(FZ_700_mz_dFz.SA*to_deg,FZ_700_mz_dFz.MZ,'o', 'Color', [0.5 1 0.5])
plot(FZ_900_mz_dFz.SA*to_deg,FZ_900_mz_dFz.MZ,'o', 'Color', [0.8 0.6 1])
plot(FZ_1120_mz_dFz.SA*to_deg,FZ_1120_mz_dFz.MZ,'o', 'Color', [0.6 1 1])
plot(SA_vec*to_deg,MZ0_fz_var_vec1,'.','LineWidth',2, 'Color', 'r')
plot(SA_vec*to_deg,MZ0_fz_var_vec2,'.','LineWidth',2, 'Color', 'b')
plot(SA_vec*to_deg,MZ0_fz_var_vec3,'.','LineWidth',2, 'Color', 'g')
plot(SA_vec*to_deg,MZ0_fz_var_vec4,'.','LineWidth',2, 'Color', 'm')
plot(SA_vec*to_deg,MZ0_fz_var_vec5,'.','LineWidth',2, 'Color', 'c')
hold off
legend({'Raw with $Fz=220N$','Raw with $Fz=440N$','Raw with $Fz=700N$','Raw with $Fz=900N$','Raw with $Fz=1120N$', '$Mz0(Fz_{220})$','$Mz0(Fz_{440})$','$Mz0(Fz_{700})$','$Mz0(Fz_{900})$','$Mz0(Fz_{1120})$'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}(Fz)$ [Nm]')

%% ---MZ0(gamma): fitting with variable camber(gamma)
% extract data with the same vertical load (Fz = 220N) 
TDataGamma_mz = FZ_220_mz;


% Extract points at constant camber and plot
GAMMA_tol_mz_dgamma = 0.05*to_rad;
idx_mz_dgamma.GAMMA_0 = 0.0*to_rad-GAMMA_tol_mz_dgamma < TDataGamma_mz.IA & TDataGamma_mz.IA < 0.0*to_rad+GAMMA_tol_mz_dgamma;
idx_mz_dgamma.GAMMA_1 = 1.0*to_rad-GAMMA_tol_mz_dgamma < TDataGamma_mz.IA & TDataGamma_mz.IA < 1.0*to_rad+GAMMA_tol_mz_dgamma;
idx_mz_dgamma.GAMMA_2 = 2.0*to_rad-GAMMA_tol_mz_dgamma < TDataGamma_mz.IA & TDataGamma_mz.IA < 2.0*to_rad+GAMMA_tol_mz_dgamma;
idx_mz_dgamma.GAMMA_3 = 3.0*to_rad-GAMMA_tol_mz_dgamma < TDataGamma_mz.IA & TDataGamma_mz.IA < 3.0*to_rad+GAMMA_tol_mz_dgamma;
idx_mz_dgamma.GAMMA_4 = 4.0*to_rad-GAMMA_tol_mz_dgamma < TDataGamma_mz.IA & TDataGamma_mz.IA < 4.0*to_rad+GAMMA_tol_mz_dgamma;

GAMMA_0__mz_dgamma  = TDataGamma_mz( idx_mz_dgamma.GAMMA_0, : );
GAMMA_1__mz_dgamma  = TDataGamma_mz( idx_mz_dgamma.GAMMA_1, : );
GAMMA_2__mz_dgamma  = TDataGamma_mz( idx_mz_dgamma.GAMMA_2, : );
GAMMA_3__mz_dgamma  = TDataGamma_mz( idx_mz_dgamma.GAMMA_3, : );
GAMMA_4__mz_dgamma  = TDataGamma_mz( idx_mz_dgamma.GAMMA_4, : );


% Fit the coeffs { qBz4, qBz5, qDz3, qDz4, qEz5, qDz8, qDz9, qHz3, qHz4 }
% Guess values for parameters to be optimised
P0__mz_dgamma = [1, 0, 0, 0, 0, 0, 0, 0, 0 ];

% Limits for parameters to be optimised
lb_mz_dgamma = [ ];
ub_mz_dgamma = [ ];

zeros_vec_mz_dgamma = zeros(size(TDataGamma_mz.IA));
ones_vec_mz_dgamma  = ones(size(TDataGamma_mz.IA));

ALPHA_vec_mz_dgamma = TDataGamma_mz.SA;
GAMMA_vec_mz_dgamma = TDataGamma_mz.IA; 
MZ_vec_mz_dgamma    = TDataGamma_mz.MZ;
FZ_vec_mz_dgamma    = TDataGamma_mz.FZ;

MZ0_varGamma_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec , GAMMA_vec_mz_dgamma, tyre_coeffs_pl.FZ0*ones(size(SA_vec)),tyre_coeffs_pl);

figure('Name','MZ0(gamma): guess', 'NumberTitle', 9 + last_fig_FY0)
plot(ALPHA_vec_mz_dgamma,TDataGamma_mz.MZ,'o')
hold on
plot(SA_vec,MZ0_varGamma_vec,'.','MarkerSize',5)
legend({'Raw data variable camber','Guess MZ0(gamma)'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}(\gamma)$ [N]')

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma_mz,fval,exitflag] = fmincon(@(P)resid_pure_Mz_varGamma(P,MZ_vec_mz_dgamma, ALPHA_vec_mz_dgamma,GAMMA_vec_mz_dgamma,tyre_coeffs_pl.FZ0, tyre_coeffs_pl),...
                               P0__mz_dgamma,[],[],[],[],lb_mz_dgamma,ub_mz_dgamma);

% Change tyre data with new optimal values       
tyre_coeffs_pl.qBz4 = P_varGamma_mz(1); 
tyre_coeffs_pl.qBz5 = P_varGamma_mz(2);
tyre_coeffs_pl.qDz3 = P_varGamma_mz(3);
tyre_coeffs_pl.qDz4 = P_varGamma_mz(4);
tyre_coeffs_pl.qEz5 = P_varGamma_mz(5);
tyre_coeffs_pl.qDz8 = P_varGamma_mz(6);
tyre_coeffs_pl.qDz9 = P_varGamma_mz(7);
tyre_coeffs_pl.qHz3 = P_varGamma_mz(8);
tyre_coeffs_pl.qHz4 = P_varGamma_mz(9);

tmp_zeros_mz_dgamma = zeros(size(SA_vec));
tmp_ones__mz_dgamma = ones(size(SA_vec));

MZ0_gamma_var_vec1 = MF96_MZ0_vec(tmp_zeros_mz_dgamma, SA_vec ,mean(GAMMA_0__mz_dgamma.IA)*tmp_ones__mz_dgamma, mean(TDataGamma.FZ)*tmp_ones__mz_dgamma,tyre_coeffs_pl);
MZ0_gamma_var_vec2 = MF96_MZ0_vec(tmp_zeros_mz_dgamma, SA_vec ,mean(GAMMA_1__mz_dgamma.IA)*tmp_ones__mz_dgamma, mean(TDataGamma.FZ)*tmp_ones__mz_dgamma,tyre_coeffs_pl);
MZ0_gamma_var_vec3 = MF96_MZ0_vec(tmp_zeros_mz_dgamma, SA_vec ,mean(GAMMA_2__mz_dgamma.IA)*tmp_ones__mz_dgamma, mean(TDataGamma.FZ)*tmp_ones__mz_dgamma,tyre_coeffs_pl);
MZ0_gamma_var_vec4 = MF96_MZ0_vec(tmp_zeros_mz_dgamma, SA_vec ,mean(GAMMA_3__mz_dgamma.IA)*tmp_ones__mz_dgamma, mean(TDataGamma.FZ)*tmp_ones__mz_dgamma,tyre_coeffs_pl);
MZ0_gamma_var_vec5 = MF96_MZ0_vec(tmp_zeros_mz_dgamma, SA_vec ,mean(GAMMA_4__mz_dgamma.IA)*tmp_ones__mz_dgamma, mean(TDataGamma.FZ)*tmp_ones__mz_dgamma,tyre_coeffs_pl);

figure('Name','MZ0(gamma): fitted with variable camber','NumberTitle', 10 + last_fig_FY0)
hold on
plot(GAMMA_0__mz_dgamma.SA*to_deg,GAMMA_0__mz_dgamma.MZ,'-','MarkerSize',5, 'Color', '#0072BD') %'MarkerEdgeColor','y',
plot(GAMMA_1__mz_dgamma.SA*to_deg,GAMMA_1__mz_dgamma.MZ,'-','MarkerSize',5, 'Color', '#D95319') %'MarkerEdgeColor','c',
plot(GAMMA_2__mz_dgamma.SA*to_deg,GAMMA_2__mz_dgamma.MZ,'-','MarkerSize',5, 'Color', '#EDB120') %'MarkerEdgeColor','m',
plot(GAMMA_3__mz_dgamma.SA*to_deg,GAMMA_3__mz_dgamma.MZ,'-','MarkerSize',5, 'Color', '#77AC30') %'MarkerEdgeColor','b',
plot(GAMMA_4__mz_dgamma.SA*to_deg,GAMMA_4__mz_dgamma.MZ,'-','MarkerSize',5, 'Color', '#4DBEEE') %'MarkerEdgeColor','r',
plot(SA_vec*to_deg,MZ0_gamma_var_vec1,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#0072BD')
plot(SA_vec*to_deg,MZ0_gamma_var_vec2,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#D95319')
plot(SA_vec*to_deg,MZ0_gamma_var_vec3,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#EDB120')
plot(SA_vec*to_deg,MZ0_gamma_var_vec4,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#77AC30')
plot(SA_vec*to_deg,MZ0_gamma_var_vec5,'-s','LineWidth',2,'MarkerSize',1, 'Color', '#4DBEEE')
legend({'Raw data with $\gamma = 0 deg $','Raw data with $\gamma = 1 deg $','Raw data with $\gamma = 2 deg $','Raw data with $\gamma = 3 deg $','Raw data with $\gamma = 4 deg $', 'Mz0($\gamma = 0 deg$)','Mz0($\gamma = 1 deg$)','Mz0($\gamma = 2 deg$)','Mz0($\gamma = 3 deg$)','Mz0($\gamma = 4 deg$)'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}$ [Nm]')


% Calculate the residuals with the optimal solution found above
res_Mz0_varGamma  = resid_pure_Mz_varGamma(P_varGamma_mz,MZ_vec_mz_dgamma, ALPHA_vec_mz_dgamma, GAMMA_vec_mz_dgamma, tyre_coeffs_pl.FZ0, tyre_coeffs_pl);


%% ---------------last figure MZ0---------------
last_fig_MZ0 = 10 + last_fig_FY0;


%% --Combined longitudinal force FX: dataset import

% Change dataset: load and select the region analysis
data_set = 'Hoosier_B1464run30'; % combined behaviour

fprintf('Loading dataset: ')

switch data_set
  case 'Hoosier_B1464run30'
      fprintf('for combined behaviour analysis.')
      load ([data_set_path, data_set]);

  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion (at the higher pressure)
cut_start_comb = 19030;
cut_end_comb   = 38170;


smpl_range_comb = cut_start_comb:cut_end_comb;

fprintf('\ncompleted!')

%% ---Dataset for combined behaviour: plot

figure ('Name','CombFXFY: entire raw dataset', 'NumberTitle',1+ last_fig_MZ0)
tiledlayout(6,1)

ax_list_5(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start_comb cut_start_comb],y_range,'--r')
plot([cut_end_comb cut_end_comb],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list_5(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start_comb cut_start_comb],y_range,'--r')
plot([cut_end_comb cut_end_comb],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_5(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start_comb cut_start_comb],y_range,'--r')
plot([cut_end_comb cut_end_comb],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_5(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start_comb cut_start_comb],y_range,'--r')
plot([cut_end_comb cut_end_comb],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list_5(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start_comb cut_start_comb],y_range,'--r')
plot([cut_end_comb cut_end_comb],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list_5(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start_comb cut_start_comb],y_range,'--r')
plot([cut_end_comb cut_end_comb],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list_5,'x')
%% ---Higher pressure dataset for combined behaviour: table selection and plot
% for different camber angle, vertical load and side slip angle

vec_samples_comb = 1:1:length(smpl_range_comb);

tyre_data_comb = table(); % create empty table
% store raw data in table
tyre_data_comb.SL =  SL(smpl_range_comb);
tyre_data_comb.SA = -SA(smpl_range_comb)*to_rad;    % SAE -> Adapted SAE
tyre_data_comb.FZ = -FZ(smpl_range_comb);           % SAE -> Adapted SAE
tyre_data_comb.FX =  FX(smpl_range_comb);
tyre_data_comb.FY =  FY(smpl_range_comb);   
tyre_data_comb.MZ =  MZ(smpl_range_comb);
tyre_data_comb.IA =  IA(smpl_range_comb)*to_rad;

% Extract points at constant camber angle
GAMMA_tol_comb = 0.05*to_rad;
idx_comb.GAMMA_0 = 0.0*to_rad-GAMMA_tol_comb < tyre_data_comb.IA & tyre_data_comb.IA < 0.0*to_rad+GAMMA_tol_comb;
idx_comb.GAMMA_1 = 1.0*to_rad-GAMMA_tol_comb < tyre_data_comb.IA & tyre_data_comb.IA < 1.0*to_rad+GAMMA_tol_comb;
idx_comb.GAMMA_2 = 2.0*to_rad-GAMMA_tol_comb < tyre_data_comb.IA & tyre_data_comb.IA < 2.0*to_rad+GAMMA_tol_comb;
idx_comb.GAMMA_3 = 3.0*to_rad-GAMMA_tol_comb < tyre_data_comb.IA & tyre_data_comb.IA < 3.0*to_rad+GAMMA_tol_comb;
idx_comb.GAMMA_4 = 4.0*to_rad-GAMMA_tol_comb < tyre_data_comb.IA & tyre_data_comb.IA < 4.0*to_rad+GAMMA_tol_comb;

GAMMA_0_comb  = tyre_data_comb( idx_comb.GAMMA_0, : );
GAMMA_1_comb  = tyre_data_comb( idx_comb.GAMMA_1, : );
GAMMA_2_comb  = tyre_data_comb( idx_comb.GAMMA_2, : );
GAMMA_3_comb  = tyre_data_comb( idx_comb.GAMMA_3, : );
GAMMA_4_comb  = tyre_data_comb( idx_comb.GAMMA_4, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol_comb = 100;
idx_comb.FZ_220  = 220-FZ_tol_comb < tyre_data_comb.FZ & tyre_data_comb.FZ < 220+FZ_tol_comb;
idx_comb.FZ_700  = 700-FZ_tol_comb < tyre_data_comb.FZ & tyre_data_comb.FZ < 700+FZ_tol_comb;
idx_comb.FZ_900  = 900-FZ_tol_comb < tyre_data_comb.FZ & tyre_data_comb.FZ < 900+FZ_tol_comb;
idx_comb.FZ_1120 = 1120-FZ_tol_comb < tyre_data_comb.FZ & tyre_data_comb.FZ < 1120+FZ_tol_comb;
FZ_220_comb  = tyre_data_comb( idx_comb.FZ_220, : );
FZ_700_comb  = tyre_data_comb( idx_comb.FZ_700, : );
FZ_900_comb  = tyre_data_comb( idx_comb.FZ_900, : );
FZ_1120_comb = tyre_data_comb( idx_comb.FZ_1120, : );

% Extract the points at constant side slip angle
% 0° , 3° , 6 °
SA_tol_comb = 0.5*to_rad;
idx_comb.SA_0    =  0-SA_tol_comb      < tyre_data_comb.SA & tyre_data_comb.SA < 0+SA_tol_comb;
idx_comb.SA_3 = (3*to_rad-SA_tol_comb) < tyre_data_comb.SA & tyre_data_comb.SA < 3*to_rad+SA_tol_comb;
idx_comb.SA_6 = (6*to_rad-SA_tol_comb) < tyre_data_comb.SA & tyre_data_comb.SA < 6*to_rad+SA_tol_comb;
SA_0_comb     = tyre_data_comb( idx_comb.SA_0, : );
SA_3_comb     = tyre_data_comb( idx_comb.SA_3, : );
SA_6_comb     = tyre_data_comb( idx_comb.SA_6, : );

% Plot
figure('Name','CombFXFY: higher pressure dataset with regions', 'NumberTitle', 2 + last_fig_MZ0)
tiledlayout(3,1)

ax_list_6(1) = nexttile;
plot(tyre_data_comb.IA*to_deg)
hold on
plot(vec_samples_comb(idx_comb.GAMMA_0),GAMMA_0_comb.IA*to_deg,'.');
plot(vec_samples_comb(idx_comb.GAMMA_1),GAMMA_1_comb.IA*to_deg,'.');
plot(vec_samples_comb(idx_comb.GAMMA_2),GAMMA_2_comb.IA*to_deg,'.');
plot(vec_samples_comb(idx_comb.GAMMA_3),GAMMA_3_comb.IA*to_deg,'.');
plot(vec_samples_comb(idx_comb.GAMMA_4),GAMMA_4_comb.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list_6(2) = nexttile;
plot(tyre_data_comb.FZ)
hold on
plot(vec_samples_comb(idx_comb.FZ_220),FZ_220_comb.FZ,'.');
%plot(vec_samples_comb(idx_comb.FZ_440),FZ_440_comb.FZ,'.');
plot(vec_samples_comb(idx_comb.FZ_700),FZ_700_comb.FZ,'.');
plot(vec_samples_comb(idx_comb.FZ_900),FZ_900_comb.FZ,'.');
plot(vec_samples_comb(idx_comb.FZ_1120),FZ_1120_comb.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list_6(3) = nexttile;
plot(tyre_data_comb.SA*to_deg)
hold on
plot(vec_samples_comb(idx_comb.SA_0),SA_0_comb.SA*to_deg,'.');
plot(vec_samples_comb(idx_comb.SA_3),SA_3_comb.SA*to_deg,'.');
plot(vec_samples_comb(idx_comb.SA_6),SA_6_comb.SA*to_deg,'.');
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

hold off
linkaxes(ax_list_6,'x')

%% ---FX: fitting in pure conditions (variable alpha, gamma = 0, Fz = 220N)
% choose the range with: variable side slip angle, camber angle = 0, vertical
% load = 220N (obv within the higher pressure dataset)

[TData_x_dalpha, ~] = intersect_table_data( GAMMA_0_comb, FZ_220_comb );

figure('Name','FX: dataset in pure conditions range', 'NumberTitle', 3 + last_fig_MZ0)
plot_selected_data(TData_x_dalpha);

% extract data with variable side slip

smpl_range_x_dalpha = size(TData_x_dalpha);
vec_samples_x_dalpha = 1:1:smpl_range_x_dalpha;

% Extract points at constant side slip and plot
ALPHA_tol_x_dalpha = 0.5*to_rad;
idx_x_dalpha.ALPHA_0 = 0.0*to_rad-ALPHA_tol_x_dalpha < TData_x_dalpha.SA & TData_x_dalpha.SA < 0.0*to_rad+ALPHA_tol_x_dalpha;
idx_x_dalpha.ALPHA_3 = 3.0*to_rad-ALPHA_tol_x_dalpha < TData_x_dalpha.SA & TData_x_dalpha.SA < 3.0*to_rad+ALPHA_tol_x_dalpha;
idx_x_dalpha.ALPHA_6 = 6.0*to_rad-ALPHA_tol_x_dalpha < TData_x_dalpha.SA & TData_x_dalpha.SA < 6.0*to_rad+ALPHA_tol_x_dalpha;

ALPHA_0_dalpha  = TData_x_dalpha( idx_x_dalpha.ALPHA_0, : );
ALPHA_3_dalpha  = TData_x_dalpha( idx_x_dalpha.ALPHA_3, : );
ALPHA_6_dalpha  = TData_x_dalpha( idx_x_dalpha.ALPHA_6, : );

% Plot
figure('Name','FX: Considered ranges for pure conditions', 'NumberTitle', 4 + last_fig_MZ0)
tiledlayout(4,1)
ax_list_7(1) = nexttile;
plot(TData_x_dalpha.SA*to_deg)
hold on
plot(vec_samples_x_dalpha(idx_x_dalpha.ALPHA_0),ALPHA_0_dalpha.SA*to_deg,'.');
plot(vec_samples_x_dalpha(idx_x_dalpha.ALPHA_3),ALPHA_3_dalpha.SA*to_deg,'.');
plot(vec_samples_x_dalpha(idx_x_dalpha.ALPHA_6),ALPHA_6_dalpha.SA*to_deg,'.');
title('Side slip angle')
xlabel('Samples [-]')
ylabel('[deg]')
hold off

ax_list_7(2) = nexttile;
plot(TData_x_dalpha.FX)
hold on
plot(vec_samples_x_dalpha(idx_x_dalpha.ALPHA_0),ALPHA_0_dalpha.FX,'.');
plot(vec_samples_x_dalpha(idx_x_dalpha.ALPHA_3),ALPHA_3_dalpha.FX,'.');
plot(vec_samples_x_dalpha(idx_x_dalpha.ALPHA_6),ALPHA_6_dalpha.FX,'.');
title('Longitudinal force')
xlabel('Samples [-]')
ylabel('[N]')
hold off

ax_list_7(3) = nexttile;
plot(TData_x_dalpha.FZ)
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')
linkaxes(ax_list_7,'x')

ax_list_7(4) = nexttile;
plot(TData_x_dalpha.IA)
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')
linkaxes(ax_list_7,'x')


% Fit the coeffs {rBx1, rBx2, rCx1, rHx1}

% Guess values for parameters to be optimised
%   [rBx1, rBx2, rCx1, rHx1]
P0_x_dalpha = [ -11 , -11 , 1 , 0 ];

% Limits for parameters to be optimised
lb_x_dalpha = [ -17, -16.5, -0.5, -0.015 ];
ub_x_dalpha = [ 20, 20, 10, 0.015];

zeros_vec_x_dalpha = zeros(size(TData_x_dalpha.IA));
ones_vec_x_dalpha  = ones(size(TData_x_dalpha.IA));

KAPPA_vec_x_dalpha = TData_x_dalpha.SL;
ALPHA_vec_x_dalpha = TData_x_dalpha.SA;
FX_vec_dalpha    = TData_x_dalpha.FX;
FZ_vec_dalpha    = TData_x_dalpha.FZ;


% Minimize the residual:
[P_x_dalpha,fval,exitflag] = fmincon(@(P)resid_Fx_varAlpha(P,FX_vec_dalpha, KAPPA_vec_x_dalpha, ALPHA_vec_x_dalpha, zeros_vec_x_dalpha,tyre_coeffs_pl.FZ0*ones_vec_x_dalpha, tyre_coeffs_pl),...
                               P0_x_dalpha,[],[],[],[],lb_x_dalpha,ub_x_dalpha);

% Change tyre data with new optimal values                             
    tyre_coeffs_pl.rBx1 = P_x_dalpha(1);  
    tyre_coeffs_pl.rBx2 = P_x_dalpha(2); 
    tyre_coeffs_pl.rCx1 = P_x_dalpha(3); 
    tyre_coeffs_pl.rHx1 = P_x_dalpha(4); 


tmp_zeros_dalpha = zeros(size(SL_vec));
tmp_ones_dalpha = ones(size(SL_vec));

[FX_alpha_var_vec1, Gxa_alpha_var_vec1] = MF96_FX_vec(SL_vec, mean(ALPHA_0_dalpha.SA)*tmp_ones_dalpha , tmp_zeros_dalpha, mean(ALPHA_0_dalpha.FZ)*tmp_ones_dalpha, tyre_coeffs_pl);
[FX_alpha_var_vec2, Gxa_alpha_var_vec2] = MF96_FX_vec(SL_vec, mean(ALPHA_3_dalpha.SA)*tmp_ones_dalpha , tmp_zeros_dalpha, mean(ALPHA_3_dalpha.FZ)*tmp_ones_dalpha,tyre_coeffs_pl);
[FX_alpha_var_vec3, Gxa_alpha_var_vec3] = MF96_FX_vec(SL_vec, mean(ALPHA_6_dalpha.SA)*tmp_ones_dalpha , tmp_zeros_dalpha, mean(ALPHA_6_dalpha.FZ)*tmp_ones_dalpha,tyre_coeffs_pl);

[~, Gxa_gamma_var_vec1] = MF96_FX_vec(0*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_0_dalpha.FZ)*ones(size(SA_vec)), tyre_coeffs_pl);
[~, Gxa_gamma_var_vec2] = MF96_FX_vec(0.1*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_3_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);
[~, Gxa_gamma_var_vec3] = MF96_FX_vec(0.3*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_6_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);


figure('Name','FX(kappa): fitted in pure conditions','NumberTitle', 5 + last_fig_MZ0)
hold on
plot(ALPHA_0_dalpha.SL,ALPHA_0_dalpha.FX,'.','MarkerSize',5, 'Color', '#0072BD') %'MarkerEdgeColor','y',
plot(ALPHA_3_dalpha.SL,ALPHA_3_dalpha.FX,'.','MarkerSize',5, 'Color', '#D95319') %'MarkerEdgeColor','c',
plot(ALPHA_6_dalpha.SL,ALPHA_6_dalpha.FX,'.','MarkerSize',5, 'Color', '#EDB120') %'MarkerEdgeColor','m',

plot(SL_vec,FX_alpha_var_vec1,'-s','LineWidth',1,'MarkerSize',1, 'Color', '#0072BD')
plot(SL_vec,FX_alpha_var_vec2,'-s','LineWidth',1,'MarkerSize',1, 'Color', '#D95319')
plot(SL_vec,FX_alpha_var_vec3,'-s','LineWidth',1,'MarkerSize',1, 'Color', '#EDB120')
legend({'Raw with $\alpha_0 = 0 deg $','Raw with $\alpha_3 = 3 deg $','Raw with $ \alpha_6 = 6 deg $', 'Fx($\alpha_0$)','Fy($\alpha_3$)','Fy($\alpha_6$)'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')

figure('Name','Gxa(kappa) as function of kappa','NumberTitle', 6 + last_fig_MZ0)
hold on
plot(SL_vec,Gxa_alpha_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
plot(SL_vec,Gxa_alpha_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
plot(SL_vec,Gxa_alpha_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
legend({'$ \alpha_0 = 0 deg $','$ \alpha_3 = 3 deg $','$ \alpha_6 = 6 deg $'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$G_{xa}$ [-]')

figure('Name','Gxa(alpha) as function of alpha','NumberTitle', 7 + last_fig_MZ0)
hold on
plot(SA_vec*to_deg,Gxa_gamma_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
plot(SA_vec*to_deg,Gxa_gamma_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
plot(SA_vec*to_deg,Gxa_gamma_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
legend({'$ \kappa = 0 $','$ \kappa = 0.1 $','$ \kappa = 0.3 $'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$G_{xa}$ [-]')


%% ---------------last figure FX---------------
% For figure number:
last_fig_FX = 7 + last_fig_MZ0;

%% --Combined longitudinal force FY: same dataset (combined)


%% ---FY: fitting in pure conditions (variable alpha, gamma = 0, Fz = 220N)
% choose the range with: variable side slip angle, camber angle = 0, vertical
% load = 220N (obv within the higher pressure dataset)

% Fit the coeffs {rBy1, rBy2, rBy3, rCy1, rHy1, rVy1, rVy4, rVy5, rVy6}

% Guess values for parameters to be optimised
%   [ rBy1, rBy2, rBy3, rCy1, rHy1, rVy1, rVy4, rVy5, rVy6 ]
P0_y_dalpha = [ 14 , 13, -0.5 , 0.98 , 0.03 , -0.23 , 3.8 , -0.1 , 28.4  ];

% Limits for parameters to be optimised
lb_y_dalpha = [ ];
ub_y_dalpha = [ ];

zeros_vec_y_dalpha = zeros(size(TData_x_dalpha.IA));
ones_vec_y_dalpha  = ones(size(TData_x_dalpha.IA));

KAPPA_vec_y_dalpha = TData_x_dalpha.SL;
ALPHA_vec_y_dalpha = TData_x_dalpha.SA;
FY_vec_dalpha    = TData_x_dalpha.FY;


% Minimize the residual:
[P_y_dalpha,fval,exitflag] = fmincon(@(P)resid_Fy_varAlpha(P,FY_vec_dalpha, KAPPA_vec_y_dalpha, ALPHA_vec_y_dalpha, zeros_vec_y_dalpha,tyre_coeffs_pl.FZ0*ones_vec_y_dalpha, tyre_coeffs_pl),...
                               P0_y_dalpha,[],[],[],[],lb_y_dalpha,ub_y_dalpha);

% Change tyre data with new optimal values                             
    tyre_coeffs_pl.rBy1 = P_y_dalpha(1);  
    tyre_coeffs_pl.rBy2 = P_y_dalpha(2); 
    tyre_coeffs_pl.rBy3 = P_y_dalpha(3); 
    tyre_coeffs_pl.rCy1 = P_y_dalpha(4);
    tyre_coeffs_pl.rHy1 = P_y_dalpha(5);
    tyre_coeffs_pl.rVy1 = P_y_dalpha(6);
    tyre_coeffs_pl.rVy4 = P_y_dalpha(7);
    tyre_coeffs_pl.rVy5 = P_y_dalpha(8);
    tyre_coeffs_pl.rVy6 = P_y_dalpha(9); 


tmp_zeros_dalpha = zeros(size(SL_vec));
tmp_ones_dalpha = ones(size(SL_vec));

[FY_alpha_var_vec1, Gyk_alpha_var_vec1] = MF96_FY_vec(SL_vec, mean(ALPHA_0_dalpha.SA)*tmp_ones_dalpha , tmp_zeros_dalpha, mean(ALPHA_0_dalpha.FZ)*tmp_ones_dalpha, tyre_coeffs_pl);
[FY_alpha_var_vec2, Gyk_alpha_var_vec2] = MF96_FY_vec(SL_vec, mean(ALPHA_3_dalpha.SA)*tmp_ones_dalpha , tmp_zeros_dalpha, mean(ALPHA_3_dalpha.FZ)*tmp_ones_dalpha,tyre_coeffs_pl);
[FY_alpha_var_vec3, Gyk_alpha_var_vec3] = MF96_FY_vec(SL_vec, mean(ALPHA_6_dalpha.SA)*tmp_ones_dalpha , tmp_zeros_dalpha, mean(ALPHA_6_dalpha.FZ)*tmp_ones_dalpha,tyre_coeffs_pl);

[~, Gyk_gamma_var_vec1] = MF96_FY_vec(0*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_0_dalpha.FZ)*ones(size(SA_vec)), tyre_coeffs_pl);
[~, Gyk_gamma_var_vec2] = MF96_FY_vec(0.1*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_3_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);
[~, Gyk_gamma_var_vec3] = MF96_FY_vec(0.3*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_6_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);

figure('Name','FY(kappa): fitted in pure conditions','NumberTitle', 1 + last_fig_FX)
hold on
plot(ALPHA_0_dalpha.SL,ALPHA_0_dalpha.FY,'.','MarkerSize',5, 'Color', '#0072BD') %'MarkerEdgeColor','y',
plot(ALPHA_3_dalpha.SL,ALPHA_3_dalpha.FY,'.','MarkerSize',5, 'Color', '#D95319') %'MarkerEdgeColor','c',
plot(ALPHA_6_dalpha.SL,ALPHA_6_dalpha.FY,'.','MarkerSize',5, 'Color', '#EDB120') %'MarkerEdgeColor','m',

plot(SL_vec,FY_alpha_var_vec1,'-s','LineWidth',1,'MarkerSize',1, 'Color', '#0072BD')
plot(SL_vec,FY_alpha_var_vec2,'-s','LineWidth',1,'MarkerSize',1, 'Color', '#D95319')
plot(SL_vec,FY_alpha_var_vec3,'-s','LineWidth',1,'MarkerSize',1, 'Color', '#EDB120')
legend({'$ \alpha_0 = 0 deg $','$ \alpha_3 = 3 deg $','$ \alpha_6 = 6 deg $', 'Fy($\alpha_0$)','Fy($\alpha_3$)','Fy($\alpha_6$)'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')

figure('Name','Gyk(kappa) as function of kappa','NumberTitle', 2 + last_fig_FX)
hold on
plot(SL_vec,Gyk_alpha_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
plot(SL_vec,Gyk_alpha_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
plot(SL_vec,Gyk_alpha_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
legend({'$ \alpha_0 = 0 deg $','$ \alpha_3 = 3 deg $','$ \alpha_6 = 6 deg $'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$G_{yk}$ [-]')

% To be done!
figure('Name','Gyk(alpha) as function of alpha','NumberTitle', 3 + last_fig_FX)
hold on
plot(SA_vec*to_deg,Gyk_gamma_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
plot(SA_vec*to_deg,Gyk_gamma_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
plot(SA_vec*to_deg,Gyk_gamma_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
legend({'$ \kappa = 0 $','$ \kappa = 0.1 $','$ \kappa = 0.3 $'}, 'Location','eastoutside');
xlabel('$\alpha$ [deg]')
ylabel('$G_{yk}$ [-]')


%% ---FY(Fz): fitting with variable Fz
% Consider the 4 cases of different vertical load and camber angle = 0, obv
% variable alpha

TData_y_comb_dFz = GAMMA_0_comb;

figure('Name','FY(Fz): considered dataset', 'NumberTitle', 4 + last_fig_FX)
plot_selected_data(TData_y_comb_dFz);

%Extract various load conditions
% =data with 0 camber and variable load

smpl_range_y_comb_dFz = size(TData_y_comb_dFz);
vec_samples_y_comb_dFz = 1:1:smpl_range_y_comb_dFz;

% Extract points at constant side slip and plot
ALPHA_tol_y_comb_dFz = 0.5*to_rad;
idx_y_comb_dFz.ALPHA_0 = 0.0*to_rad-ALPHA_tol_y_comb_dFz < TData_y_comb_dFz.SA & TData_y_comb_dFz.SA < 0.0*to_rad+ALPHA_tol_y_comb_dFz;
idx_y_comb_dFz.ALPHA_3 = 3.0*to_rad-ALPHA_tol_y_comb_dFz < TData_y_comb_dFz.SA & TData_y_comb_dFz.SA < 3.0*to_rad+ALPHA_tol_y_comb_dFz;
idx_y_comb_dFz.ALPHA_6 = 6.0*to_rad-ALPHA_tol_y_comb_dFz < TData_y_comb_dFz.SA & TData_y_comb_dFz.SA < 6.0*to_rad+ALPHA_tol_y_comb_dFz;

ALPHA_0_y_comb_dFz  = TData_y_comb_dFz( idx_y_comb_dFz.ALPHA_0, : );
ALPHA_3_y_comb_dFz  = TData_y_comb_dFz( idx_y_comb_dFz.ALPHA_3, : );
ALPHA_6_y_comb_dFz  = TData_y_comb_dFz( idx_y_comb_dFz.ALPHA_6, : );


FZ_tol_y_comb_dFz = 100;
idx_y_comb_dFz.FZ_220  = 220-FZ_tol_y_comb_dFz < TData_y_comb_dFz.FZ & TData_y_comb_dFz.FZ < 220+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.FZ_700  = 700-FZ_tol_y_comb_dFz < TData_y_comb_dFz.FZ & TData_y_comb_dFz.FZ < 700+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.FZ_900  = 900-FZ_tol_y_comb_dFz < TData_y_comb_dFz.FZ & TData_y_comb_dFz.FZ < 900+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.FZ_1120 = 1120-FZ_tol_y_comb_dFz < TData_y_comb_dFz.FZ & TData_y_comb_dFz.FZ < 1120+FZ_tol_y_comb_dFz;
FZ_220_y_comb_dFz  = TData_y_comb_dFz( idx_y_comb_dFz.FZ_220, : );
FZ_700_y_comb_dFz  = TData_y_comb_dFz( idx_y_comb_dFz.FZ_700, : );
FZ_900_y_comb_dFz  = TData_y_comb_dFz( idx_y_comb_dFz.FZ_900, : );
FZ_1120_y_comb_dFz = TData_y_comb_dFz( idx_y_comb_dFz.FZ_1120, : );

% Ranges for plot (Skip)
idx_y_comb_dFz.alpha0_FZ_220  = 220-FZ_tol_y_comb_dFz < ALPHA_0_y_comb_dFz.FZ & ALPHA_0_y_comb_dFz.FZ < 220+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha0_FZ_700  = 700-FZ_tol_y_comb_dFz < ALPHA_0_y_comb_dFz.FZ & ALPHA_0_y_comb_dFz.FZ < 700+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha0_FZ_900  = 900-FZ_tol_y_comb_dFz < ALPHA_0_y_comb_dFz.FZ & ALPHA_0_y_comb_dFz.FZ < 900+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha0_FZ_1120 = 1120-FZ_tol_y_comb_dFz < ALPHA_0_y_comb_dFz.FZ & ALPHA_0_y_comb_dFz.FZ < 1120+FZ_tol_y_comb_dFz;
alpha0_FZ_220_y_comb_dFz  = ALPHA_0_y_comb_dFz( idx_y_comb_dFz.alpha0_FZ_220, : );
alpha0_FZ_700_y_comb_dFz  = ALPHA_0_y_comb_dFz( idx_y_comb_dFz.alpha0_FZ_700, : );
alpha0_FZ_900_y_comb_dFz  = ALPHA_0_y_comb_dFz( idx_y_comb_dFz.alpha0_FZ_900, : );
alpha0_FZ_1120_y_comb_dFz = ALPHA_0_y_comb_dFz( idx_y_comb_dFz.alpha0_FZ_1120, : );

idx_y_comb_dFz.alpha3_FZ_220  = 220-FZ_tol_y_comb_dFz < ALPHA_3_y_comb_dFz.FZ & ALPHA_3_y_comb_dFz.FZ < 220+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha3_FZ_700  = 700-FZ_tol_y_comb_dFz < ALPHA_3_y_comb_dFz.FZ & ALPHA_3_y_comb_dFz.FZ < 700+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha3_FZ_900  = 900-FZ_tol_y_comb_dFz < ALPHA_3_y_comb_dFz.FZ & ALPHA_3_y_comb_dFz.FZ < 900+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha3_FZ_1120 = 1120-FZ_tol_y_comb_dFz < ALPHA_3_y_comb_dFz.FZ & ALPHA_3_y_comb_dFz.FZ < 1120+FZ_tol_y_comb_dFz;
alpha3_FZ_220_y_comb_dFz  = ALPHA_3_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_220, : );
alpha3_FZ_700_y_comb_dFz  = ALPHA_3_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_700, : );
alpha3_FZ_900_y_comb_dFz  = ALPHA_3_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_900, : );
alpha3_FZ_1120_y_comb_dFz = ALPHA_3_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_1120, : );

idx_y_comb_dFz.alpha6_FZ_220  = 220-FZ_tol_y_comb_dFz < ALPHA_6_y_comb_dFz.FZ & ALPHA_6_y_comb_dFz.FZ < 220+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha6_FZ_700  = 700-FZ_tol_y_comb_dFz < ALPHA_6_y_comb_dFz.FZ & ALPHA_6_y_comb_dFz.FZ < 700+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha6_FZ_900  = 900-FZ_tol_y_comb_dFz < ALPHA_6_y_comb_dFz.FZ & ALPHA_6_y_comb_dFz.FZ < 900+FZ_tol_y_comb_dFz;
idx_y_comb_dFz.alpha6_FZ_1120 = 1120-FZ_tol_y_comb_dFz < ALPHA_6_y_comb_dFz.FZ & ALPHA_6_y_comb_dFz.FZ < 1120+FZ_tol_y_comb_dFz;
alpha6_FZ_220_y_comb_dFz  = ALPHA_6_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_220, : );
alpha6_FZ_700_y_comb_dFz  = ALPHA_6_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_700, : );
alpha6_FZ_900_y_comb_dFz  = ALPHA_6_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_900, : );
alpha6_FZ_1120_y_comb_dFz = ALPHA_6_y_comb_dFz( idx_y_comb_dFz.alpha3_FZ_1120, : );

% Plot
figure('Name','FY(Fz): dataset with regions', 'NumberTitle', 5 + last_fig_FX)
tiledlayout(3,1)
ax_list_8(1) = nexttile;
plot(TData_y_comb_dFz.SA*to_deg)
hold on
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.ALPHA_0),ALPHA_0_y_comb_dFz.SA*to_deg,'.');
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.ALPHA_3),ALPHA_3_y_comb_dFz.SA*to_deg,'.');
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.ALPHA_6),ALPHA_6_y_comb_dFz.SA*to_deg,'.');
title('Side slip angle')
xlabel('Samples [-]')
ylabel('[deg]')
hold off

ax_list_8(2) = nexttile;
plot(TData_y_comb_dFz.FY)
hold on
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.ALPHA_0),ALPHA_0_y_comb_dFz.FY,'.');
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.ALPHA_3),ALPHA_3_y_comb_dFz.FY,'.');
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.ALPHA_6),ALPHA_6_y_comb_dFz.FY,'.');
title('Lateral force')
xlabel('Samples [-]')
ylabel('[N]')
hold off

ax_list_8(3) = nexttile;
plot(TData_y_comb_dFz.FZ)
hold on
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.FZ_220),FZ_220_y_comb_dFz.FZ,'.');
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.FZ_700),FZ_700_y_comb_dFz.FZ,'.');
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.FZ_900),FZ_900_y_comb_dFz.FZ,'.');
plot(vec_samples_y_comb_dFz(idx_y_comb_dFz.FZ_1120),FZ_1120_y_comb_dFz.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')
hold off
linkaxes(ax_list_8,'x')

% Fit the coeffs {rVy2}

% Guess values for parameters to be optimised
%   [rVy2]
P0_y_comb_dFz = [0];

% Limits for parameters to be optimised
lb_y_comb_dFz = [ ];
ub_y_comb_dFz = [ ];

zeros_vec_y_comb_dFz = zeros(size(TData_y_comb_dFz.IA));
ones_vec_y_comb_dFz  = ones(size(TData_y_comb_dFz.IA));

KAPPA_vec_y_comb_dFz = TData_y_comb_dFz.SL;
ALPHA_vec_y_comb_dFz = TData_y_comb_dFz.SA;
FY_vec_comb_dFz    = TData_y_comb_dFz.FY;
FZ_vec_comb_dFz    = TData_y_comb_dFz.FZ;


% LSM_pure_Fx returns the residual, so minimize the residual varying alpha:
[P_y_comb_dFz,fval,exitflag] = fmincon(@(P)resid_Fy_varFz(P,FY_vec_comb_dFz, KAPPA_vec_y_comb_dFz, ALPHA_vec_y_comb_dFz, zeros_vec_y_comb_dFz,FZ_vec_comb_dFz, tyre_coeffs_pl),...
                               P0_y_comb_dFz,[],[],[],[],lb_y_comb_dFz,ub_y_comb_dFz);

% Change tyre data with new optimal values                             
    tyre_coeffs_pl.rVy2 = P_y_comb_dFz(1);  
 

[FY_comb_dFz_vec,~] = MF96_FY_vec(KAPPA_vec_y_comb_dFz, ALPHA_vec_y_comb_dFz, zeros_vec_y_comb_dFz, tyre_coeffs_pl.FZ0*ones_vec_y_comb_dFz, tyre_coeffs_pl);

% figure('Name','Fx vs Alpha', 'NumberTitle', 6 + last_fig_FY0)
% plot(KAPPA_vec_y_comb_dFz,TData_y_comb_dFz.FY,'o')
% hold on
% plot(KAPPA_vec_y_comb_dFz,FY_comb_dFz_vec,'-')
% xlabel('$\kappa$ [-]')
% ylabel('$F_{y} (Fz)$ [N]')

tmp_zeros_comb_dFz = zeros(size(SL_vec));
tmp_ones_comb_dFz = ones(size(SL_vec));

[FY_dFz_var_vec1, Gxa_dFz_var_vec1] = MF96_FY_vec(SL_vec, mean(ALPHA_0_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(ALPHA_0_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);
[FY_dFz_var_vec2, Gxa_dFz_var_vec2] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(ALPHA_3_y_comb_dFz.FZ)*tmp_ones_comb_dFz,tyre_coeffs_pl);
[FY_dFz_var_vec3, Gxa_dFz_var_vec3] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(ALPHA_6_y_comb_dFz.FZ)*tmp_ones_comb_dFz,tyre_coeffs_pl);

% [~, Gxa_gamma_var_vec1] = MF96_FX_vec(0*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_0_dalpha.FZ)*ones(size(SA_vec)), tyre_coeffs_pl);
% [~, Gxa_gamma_var_vec2] = MF96_FX_vec(0.1*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_3_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);
% [~, Gxa_gamma_var_vec3] = MF96_FX_vec(0.2*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_6_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);


figure('Name','FY(Fz) for all side slip angles','NumberTitle', 7 + last_fig_FX)
hold on
plot(ALPHA_0_y_comb_dFz.SL,ALPHA_0_y_comb_dFz.FY,'.','MarkerSize',5) %'MarkerEdgeColor','y',
plot(ALPHA_3_y_comb_dFz.SL,ALPHA_3_y_comb_dFz.FY,'.','MarkerSize',5) %'MarkerEdgeColor','c',
plot(ALPHA_6_y_comb_dFz.SL,ALPHA_6_y_comb_dFz.FY,'.','MarkerSize',5) %'MarkerEdgeColor','m',
plot(SL_vec,FY_dFz_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
plot(SL_vec,FY_dFz_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
plot(SL_vec,FY_dFz_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
legend({'$ \alpha_0 = 0 deg $','$ \alpha_3 = 3 deg $','$ \alpha_6 = 6 deg $', 'Fx($\alpha_0$)','Fy($\alpha_3$)','Fy($\alpha_6$)'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{y}$ [N]')


[FY_dFz220_alpha3_vec, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_220_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);
[FY_dFz700_alpha3_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_700_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);
[FY_dFz900_alpha3_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_900_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);
[FY_dFz1120_alpha3_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_1120_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);

figure('Name','FY(Fz) with alpha = 3 deg','NumberTitle', 8 + last_fig_FX)
hold on
plot(alpha3_FZ_220_y_comb_dFz.SL,alpha3_FZ_220_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[1 0.75 0.79]) %'MarkerEdgeColor','m',
plot(alpha3_FZ_700_y_comb_dFz.SL,alpha3_FZ_700_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[1 0.85 0.64]) %'MarkerEdgeColor','m',
plot(alpha3_FZ_900_y_comb_dFz.SL,alpha3_FZ_900_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[0.79 0.88 1]) %'MarkerEdgeColor','m',
plot(alpha3_FZ_1120_y_comb_dFz.SL,alpha3_FZ_1120_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[0.85 0.75 1]) %'MarkerEdgeColor','m',
plot(SL_vec,FY_dFz220_alpha3_vec,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 0.75 0.79])
plot(SL_vec,FY_dFz700_alpha3_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 0.7 0.34])
plot(SL_vec,FY_dFz900_alpha3_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[0.64 0.75 1])
plot(SL_vec,FY_dFz1120_alpha3_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[0.7 0.5 1])
legend({'Raw with $Fz=220N$,$\alpha=3 deg$','Raw with $Fz=700N$,$\alpha=3 deg$','Raw with $Fz=900N$,$\alpha=3 deg$','Raw with $Fz=1120N$,$\alpha=3 deg$','Fy(Fz=220N), fitted','Fy(Fz=700N), fitted', 'Fy(Fz=900N), fitted','Fy(Fz=1120N), fitted'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{y}(Fz)$ [N]')

[FY_dFz220_alpha6_vec, Gxa_dFz220_var_vec1] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_220_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);
[FY_dFz700_alpha6_vec1, Gxa_dFz700_var_vec1] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_700_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);
[FY_dFz900_alpha6_vec1, Gxa_dFz900_var_vec1] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_900_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);
[FY_dFz1120_alpha6_vec1, Gxa_dFz1120_var_vec1] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dFz.SA)*tmp_ones_comb_dFz , tmp_zeros_comb_dFz, mean(FZ_1120_y_comb_dFz.FZ)*tmp_ones_comb_dFz, tyre_coeffs_pl);

figure('Name','FY(Fz) with alpha = 6 deg','NumberTitle', 9 + last_fig_FX)
hold on
plot(alpha6_FZ_220_y_comb_dFz.SL,alpha6_FZ_220_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[1 0.75 0.79]) %'MarkerEdgeColor','m',
plot(alpha6_FZ_700_y_comb_dFz.SL,alpha6_FZ_700_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[1 0.85 0.64]) %'MarkerEdgeColor','m',
plot(alpha6_FZ_900_y_comb_dFz.SL,alpha6_FZ_900_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[0.79 0.88 1]) %'MarkerEdgeColor','m',
plot(alpha6_FZ_1120_y_comb_dFz.SL,alpha6_FZ_1120_y_comb_dFz.FY,'.','MarkerSize',5,'Color',[0.85 0.75 1]) %'MarkerEdgeColor','m',
plot(SL_vec,FY_dFz220_alpha6_vec,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 0.75 0.79])
plot(SL_vec,FY_dFz700_alpha6_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 0.7 0.34])
plot(SL_vec,FY_dFz900_alpha6_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[0.64 0.75 1])
plot(SL_vec,FY_dFz1120_alpha6_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[0.7 0.5 1])
legend({'Raw with $Fz=220N$,$\alpha=6 deg$','Raw with $Fz=700N$,$\alpha=6 deg$','Raw with $Fz=900N$,$\alpha=6 deg$','Raw with $Fz=1120N$,$\alpha=6 deg$','Fy(Fz=220N), fitted','Fy(Fz=700N), fitted', 'Fy(Fz=900N), fitted','Fy(Fz=1120N), fitted'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{y}(Fz)$ [N]')


% figure('Name','Gxa coeffs as function of kappa','NumberTitle', 8 + last_fig_FY0)
% hold on
% plot(SL_vec,Gxa_dFz_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
% plot(SL_vec,Gxa_dFz_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
% plot(SL_vec,Gxa_dFz_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
% legend({'$ \alpha_0 = 0 deg $','$ \alpha_3 = 3 deg $','$ \alpha_6 = 6 deg $'}, 'Location','eastoutside');
% xlabel('$\kappa$ [-]')
% ylabel('$G_{xa}$ [-]')
% 
% % To be done!
% % figure('Name','Gxa coeffs as function of alpha','NumberTitle', 9 + last_fig_FY0)
% % hold on
% % plot(SA_vec,Gxa_gamma_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
% % plot(SA_vec,Gxa_gamma_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
% % plot(SA_vec,Gxa_gamma_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
% % legend({'$ \kappa = 0 $','$ \kappa = 0.1 $','$ \kappa = 0.3 $'}, 'Location','eastoutside');
% % xlabel('$\alpha$ [rad]')
% % ylabel('$G_{xa}$ [-]')

%% ---FY(gamma): fitting with variable camber (gamma)
% evaluate the differences at the same nominal load Fz = 220N

TData_y_comb_dgamma = FZ_220_comb;

smpl_range_y_comb_dgamma = size(TData_y_comb_dgamma);
vec_samples_y_comb_dgamma = 1:1:smpl_range_y_comb_dgamma;

% Extract points at constant side slip and plot
ALPHA_tol_y_comb_dgamma = 0.5*to_rad;
idx_y_comb_dgamma.ALPHA_0 = 0.0*to_rad-ALPHA_tol_y_comb_dgamma < TData_y_comb_dgamma.SA & TData_y_comb_dgamma.SA < 0.0*to_rad+ALPHA_tol_y_comb_dgamma;
idx_y_comb_dgamma.ALPHA_3 = 3.0*to_rad-ALPHA_tol_y_comb_dgamma < TData_y_comb_dgamma.SA & TData_y_comb_dgamma.SA < 3.0*to_rad+ALPHA_tol_y_comb_dgamma;
idx_y_comb_dgamma.ALPHA_6 = 6.0*to_rad-ALPHA_tol_y_comb_dgamma < TData_y_comb_dgamma.SA & TData_y_comb_dgamma.SA < 6.0*to_rad+ALPHA_tol_y_comb_dgamma;

ALPHA_0_y_comb_dgamma  = TData_y_comb_dgamma( idx_y_comb_dgamma.ALPHA_0, : );
ALPHA_3_y_comb_dgamma  = TData_y_comb_dgamma( idx_y_comb_dgamma.ALPHA_3, : );
ALPHA_6_y_comb_dgamma  = TData_y_comb_dgamma( idx_y_comb_dgamma.ALPHA_6, : );

% Extract points at constant camber and plot
GAMMA_tol_y_comb_dgamma = 0.05*to_rad;
idx_y_comb_dgamma.GAMMA_0 = 0.0*to_rad-GAMMA_tol_y_comb_dgamma < TData_y_comb_dgamma.IA & TData_y_comb_dgamma.IA < 0.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.GAMMA_1 = 2.0*to_rad-GAMMA_tol_y_comb_dgamma < TData_y_comb_dgamma.IA & TData_y_comb_dgamma.IA < 2.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.GAMMA_2 = 4.0*to_rad-GAMMA_tol_y_comb_dgamma < TData_y_comb_dgamma.IA & TData_y_comb_dgamma.IA < 4.0*to_rad+GAMMA_tol_y_comb_dgamma;

GAMMA_0_y_comb_dgamma  = TData_y_comb_dgamma( idx_y_comb_dgamma.GAMMA_0, : );
GAMMA_1_y_comb_dgamma  = TData_y_comb_dgamma( idx_y_comb_dgamma.GAMMA_1, : );
GAMMA_2_y_comb_dgamma  = TData_y_comb_dgamma( idx_y_comb_dgamma.GAMMA_2, : );

% Ranges for plot (Skip)
idx_y_comb_dgamma.alpha0_gamma0  = 0.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_0_y_comb_dgamma.IA & ALPHA_0_y_comb_dgamma.IA < 0.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.alpha0_gamma2  = 2.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_0_y_comb_dgamma.IA & ALPHA_0_y_comb_dgamma.IA < 2.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.alpha0_gamma4  = 4.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_0_y_comb_dgamma.IA & ALPHA_0_y_comb_dgamma.IA < 4.0*to_rad+GAMMA_tol_y_comb_dgamma;
alpha0_gamma0_y_comb_dgamma  = ALPHA_0_y_comb_dgamma( idx_y_comb_dgamma.alpha0_gamma0, : );
alpha0_gamma2_y_comb_dgamma  = ALPHA_0_y_comb_dgamma( idx_y_comb_dgamma.alpha0_gamma2, : );
alpha0_gamma4_y_comb_dgamma  = ALPHA_0_y_comb_dgamma( idx_y_comb_dgamma.alpha0_gamma4, : );

idx_y_comb_dgamma.alpha3_gamma0  = 0.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_3_y_comb_dgamma.IA & ALPHA_3_y_comb_dgamma.IA < 0.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.alpha3_gamma2  = 2.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_3_y_comb_dgamma.IA & ALPHA_3_y_comb_dgamma.IA < 2.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.alpha3_gamma4  = 4.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_3_y_comb_dgamma.IA & ALPHA_3_y_comb_dgamma.IA < 4.0*to_rad+GAMMA_tol_y_comb_dgamma;
alpha3_gamma0_y_comb_dgamma  = ALPHA_3_y_comb_dgamma( idx_y_comb_dgamma.alpha3_gamma0, : );
alpha3_gamma2_y_comb_dgamma  = ALPHA_3_y_comb_dgamma( idx_y_comb_dgamma.alpha3_gamma2, : );
alpha3_gamma4_y_comb_dgamma  = ALPHA_3_y_comb_dgamma( idx_y_comb_dgamma.alpha3_gamma4, : );

idx_y_comb_dgamma.alpha6_gamma0  = 0.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_6_y_comb_dgamma.IA & ALPHA_6_y_comb_dgamma.IA < 0.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.alpha6_gamma2  = 2.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_6_y_comb_dgamma.IA & ALPHA_6_y_comb_dgamma.IA < 2.0*to_rad+GAMMA_tol_y_comb_dgamma;
idx_y_comb_dgamma.alpha6_gamma4  = 4.0*to_rad-GAMMA_tol_y_comb_dgamma < ALPHA_6_y_comb_dgamma.IA & ALPHA_6_y_comb_dgamma.IA < 4.0*to_rad+GAMMA_tol_y_comb_dgamma;
alpha6_gamma0_y_comb_dgamma  = ALPHA_6_y_comb_dgamma( idx_y_comb_dgamma.alpha6_gamma0, : );
alpha6_gamma2_y_comb_dgamma  = ALPHA_6_y_comb_dgamma( idx_y_comb_dgamma.alpha6_gamma2, : );
alpha6_gamma4_y_comb_dgamma  = ALPHA_6_y_comb_dgamma( idx_y_comb_dgamma.alpha6_gamma4, : );

% Plot
figure('Name','FY(gamma): dataset with regions', 'NumberTitle', 10 + last_fig_FX)
tiledlayout(3,1)
ax_list_9(1) = nexttile;
plot(TData_y_comb_dgamma.IA*to_deg)
hold on
plot(vec_samples_y_comb_dgamma(idx_y_comb_dgamma.GAMMA_0),GAMMA_0_y_comb_dgamma.IA*to_deg,'.');
plot(vec_samples_y_comb_dgamma(idx_y_comb_dgamma.GAMMA_1),GAMMA_1_y_comb_dgamma.IA*to_deg,'.');
plot(vec_samples_y_comb_dgamma(idx_y_comb_dgamma.GAMMA_2),GAMMA_2_y_comb_dgamma.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')
hold off

ax_list_9(2) = nexttile;
plot(TData_y_comb_dgamma.FY)
hold on
plot(vec_samples_y_comb_dgamma(idx_y_comb_dgamma.GAMMA_0),GAMMA_0_y_comb_dgamma.FY,'.');
plot(vec_samples_y_comb_dgamma(idx_y_comb_dgamma.GAMMA_1),GAMMA_1_y_comb_dgamma.FY,'.');
plot(vec_samples_y_comb_dgamma(idx_y_comb_dgamma.GAMMA_2),GAMMA_2_y_comb_dgamma.FY,'.');
title('Lateral force')
xlabel('Samples [-]')
ylabel('[N]')
hold off

ax_list_9(3) = nexttile;
plot(TData_y_comb_dgamma.FZ)
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')
linkaxes(ax_list_9,'x')

% Fit the coeffs {rVy3}

% Guess values for parameters to be optimised
%   [rVy3]
P0_y_comb_dgamma = [0];

% Limits for parameters to be optimised
lb_y_comb_dgamma = [ ];
ub_y_comb_dgamma = [ ];

zeros_vec_y_comb_dgamma = zeros(size(TData_y_comb_dgamma.IA));
ones_vec_y_comb_dgamma  = ones(size(TData_y_comb_dgamma.IA));

GAMMA_vec_y_comb_dgamma = TData_y_comb_dgamma.IA;
ALPHA_vec_y_comb_dgamma = TData_y_comb_dgamma.SA; 
KAPPA_vec_y_comb_dgamma = TData_y_comb_dgamma.SL; 
FY_vec_y_comb_dgamma    = TData_y_comb_dgamma.FY;
FZ_vec_y_comb_dgamma    = TData_y_comb_dgamma.FZ;


% Minimization of residuals
[P_y_comb_varGamma,fval,exitflag] = fmincon(@(P)resid_Fy_varGamma(P,FY_vec_y_comb_dgamma, KAPPA_vec_y_comb_dgamma, ALPHA_vec_y_comb_dgamma, GAMMA_vec_y_comb_dgamma,mean(FZ_vec_y_comb_dgamma)*ones_vec_y_comb_dgamma, tyre_coeffs_pl),...
                               P0_y_comb_dgamma,[],[],[],[],lb_y_comb_dgamma,ub_y_comb_dgamma);

% Change tyre data with new optimal values                             
tyre_coeffs_pl.rVy3 = P_y_comb_varGamma(1);  

tmp_zeros_comb_dgamma = zeros(size(SL_vec));
tmp_ones_comb_dgamma = ones(size(SL_vec));

% [FY_dgamma0_var_vec1, Gxa_dgamma_var_vec1] = MF96_FY_vec(SL_vec, mean(ALPHA_0_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_0_y_comb_dgamma.SA)*tmp_ones_comb_dgamma, mean(ALPHA_0_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma, tyre_coeffs_pl);
% [FY_dgamma_var_vec2, Gxa_dgamma_var_vec2] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_0_y_comb_dgamma.SA)*tmp_ones_comb_dgamma, mean(ALPHA_0_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma,tyre_coeffs_pl);
% [FY_dgamma_var_vec3, Gxa_dgamma_var_vec3] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_0_y_comb_dgamma.SA)*tmp_ones_comb_dgamma, mean(ALPHA_0_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma,tyre_coeffs_pl);
% 
% % [~, Gxa_gamma_var_vec1] = MF96_FX_vec(0*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_0_dalpha.FZ)*ones(size(SA_vec)), tyre_coeffs_pl);
% % [~, Gxa_gamma_var_vec2] = MF96_FX_vec(0.1*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_3_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);
% % [~, Gxa_gamma_var_vec3] = MF96_FX_vec(0.2*ones(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), mean(ALPHA_6_dalpha.FZ)*ones(size(SA_vec)),tyre_coeffs_pl);
% 
% figure('Name','FY(Fz) for all side slip angles','NumberTitle', 22 + last_fig_FY0)
% hold on
% plot(ALPHA_0_y_comb_dgamma.SL,ALPHA_0_y_comb_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','y',
% plot(ALPHA_3_y_comb_dgamma.SL,ALPHA_3_y_comb_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','c',
% plot(ALPHA_6_y_comb_dgamma.SL,ALPHA_6_y_comb_dgamma.FY,'.','MarkerSize',5) %'MarkerEdgeColor','m',
% plot(SL_vec,FY_dgamma0_var_vec1,'-s','LineWidth',1,'MarkerSize',1)
% plot(SL_vec,FY_dgamma_var_vec2,'-s','LineWidth',1,'MarkerSize',1)
% plot(SL_vec,FY_dgamma_var_vec3,'-s','LineWidth',1,'MarkerSize',1)
% legend({'$ \alpha_0 = 0 deg $','$ \alpha_3 = 3 deg $','$ \alpha_6 = 6 deg $', 'Fy($\alpha_0$)','Fy($\alpha_3$)','Fy($\alpha_6$)'}, 'Location','eastoutside');
% xlabel('$\kappa$ [-]')
% ylabel('$F_{y}$ [N]')

[FY_dgamma0_alpha3_var_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_0_y_comb_dgamma.IA)*tmp_ones_comb_dgamma, mean(ALPHA_3_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma, tyre_coeffs_pl);
[FY_dgamma1_alpha3_var_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_1_y_comb_dgamma.IA)*tmp_ones_comb_dgamma, mean(ALPHA_3_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma, tyre_coeffs_pl);
[FY_dgamma2_alpha3_var_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_3_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_2_y_comb_dgamma.IA)*tmp_ones_comb_dgamma, mean(ALPHA_3_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma, tyre_coeffs_pl);

figure('Name','FY(gamma) with alpha = 3deg','NumberTitle', 11 + last_fig_FX)
hold on
plot(alpha3_gamma0_y_comb_dgamma.SL,alpha3_gamma0_y_comb_dgamma.FY,'.','MarkerSize',5,'Color',[1 0.75 0.79])
plot(alpha3_gamma2_y_comb_dgamma.SL,alpha3_gamma2_y_comb_dgamma.FY,'.','MarkerSize',5,'Color',[1 1 0.4])
plot(alpha3_gamma4_y_comb_dgamma.SL,alpha3_gamma4_y_comb_dgamma.FY,'.','MarkerSize',5,'Color',[0.79 0.88 1])
plot(SL_vec,FY_dgamma0_alpha3_var_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 0.5 0.64])
plot(SL_vec,FY_dgamma1_alpha3_var_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 1 0])
plot(SL_vec,FY_dgamma2_alpha3_var_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[0.64 0.75 1])
legend({'Raw with $\gamma=0deg$,$\alpha=3 deg$','Raw with $\gamma=2deg$,$\alpha=3 deg$','Raw with $\gamma=2deg$,$\alpha=3 deg$','Fy($\gamma=0deg$), fitted','Fy($\gamma=2deg$), fitted', 'Fy($\gamma=4deg$), fitted'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{y}(\gamma)$ [N]')

[FY_dgamma0_alpha6_var_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_0_y_comb_dgamma.IA)*tmp_ones_comb_dgamma, mean(ALPHA_6_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma, tyre_coeffs_pl);
[FY_dgamma1_alpha6_var_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_1_y_comb_dgamma.IA)*tmp_ones_comb_dgamma, mean(ALPHA_6_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma, tyre_coeffs_pl);
[FY_dgamma2_alpha6_var_vec1, ~] = MF96_FY_vec(SL_vec, mean(ALPHA_6_y_comb_dgamma.SA)*tmp_ones_comb_dgamma , mean(GAMMA_2_y_comb_dgamma.IA)*tmp_ones_comb_dgamma, mean(ALPHA_6_y_comb_dgamma.FZ)*tmp_ones_comb_dgamma, tyre_coeffs_pl);

figure('Name','FY(gamma) with alpha = 6deg','NumberTitle', 12 + last_fig_FX)
hold on
plot(alpha6_gamma0_y_comb_dgamma.SL,alpha6_gamma0_y_comb_dgamma.FY,'.','MarkerSize',5,'Color',[1 0.75 0.79])
plot(alpha6_gamma2_y_comb_dgamma.SL,alpha6_gamma2_y_comb_dgamma.FY,'.','MarkerSize',5,'Color',[1 1 0.4])
plot(alpha6_gamma4_y_comb_dgamma.SL,alpha6_gamma4_y_comb_dgamma.FY,'.','MarkerSize',5,'Color',[0.79 0.88 1])
plot(SL_vec,FY_dgamma0_alpha6_var_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 0.5 0.64])
plot(SL_vec,FY_dgamma1_alpha6_var_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[1 1 0])
plot(SL_vec,FY_dgamma2_alpha6_var_vec1,'-s','LineWidth',1,'MarkerSize',1,'Color',[0.64 0.75 1])
legend({'Raw with $\gamma=0deg$,$\alpha=6 deg$','Raw with $\gamma=2deg$,$\alpha=6 deg$','Raw with $\gamma=2deg$,$\alpha=6 deg$','Fy($\gamma=0deg$), fitted','Fy($\gamma=2deg$), fitted', 'Fy($\gamma=4deg$), fitted'}, 'Location','eastoutside');
xlabel('$\kappa$ [-]')
ylabel('$F_{y}(\gamma)$ [N]')

%% -Save tyre data structure to mat file
save('tyre_coeffs_team6.mat','tyre_coeffs_pl');

