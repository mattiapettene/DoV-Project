clc
clear all
close all

format long

%% Points of the ROLL ANGLE - CAMBER GAIN curve

% Load the file. The first column contains roll angle data (in °),
% while the second column contains camber gain data (in °)
load Camber_Model.txt;

roll_angle = Camber_Model(:,1); 
camber_gain = Camber_Model(:,2);

figure(1)
plot(roll_angle,camber_gain,'ro');
grid on

%% Fitting of the model with 1st order polynomial

F = @(x,roll_angle)(x(1)*roll_angle); 
x0 = [1];

options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',50000,'MaxIterations',50000,'FunctionTolerance',1e-9,...
    'StepTolerance',1e-9,'Display','final-detailed');
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,roll_angle,camber_gain,[],[],options);

fitted_model = x(1)*camber_gain;

figure(2)
plot(roll_angle,camber_gain,'ro');
grid on
hold on
plot(roll_angle,fitted_model,'-g','LineWidth',1.5)
hold off
grid on
set(gca,'fontsize',14)  
set(gcf,'units','points','position',[50,50,460,320])
xlim([-3 3])
ylim([-2.5 3])
legend('Real','Fitted','Location','southeast')  
title('Camber gain model','interpreter','latex','FontSize',18)
xlabel('Roll angle $\phi$ (deg)','interpreter','latex','FontSize',18)
ylabel('Camber gain (deg)','interpreter','latex','FontSize',18)
set(findall(gcf,'type','line'),'linewidth',2) 






