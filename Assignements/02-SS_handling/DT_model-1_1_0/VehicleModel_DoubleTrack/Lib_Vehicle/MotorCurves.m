% clc; clear all; 
close all;

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Set defualt figure window as 'docked'
set(0,'DefaultFigureWindowStyle','docked');

%% This script displays the torque delivered to the wheel at different operating conditions

params = getVehicleDataStruct;
% Load the transmission ratio and the efficiency of the gearbox (Chimera EVO FSAE vehicle)
tau_red = params.transmission.tau_red; 
eff_red = params.transmission.eff_red;

maxTorque = params.motor.maxTorque;     % [Nm] max motor torque

maxRotSpeed = params.motor.maxRotSpeed; % [rpm] Max rotational speed of the motor 
omega_mot_max = maxRotSpeed*pi/30;

omega__wheel = 0:0.5:omega_mot_max/tau_red;

% Case 1
Tm_req_1 = maxTorque;  % 80 Nm
Tw_1 = zeros(1,length(omega__wheel));
for i=1:length(omega__wheel)
    Tw_1(i) = MotorModel(Tm_req_1,omega__wheel(i),tau_red,eff_red);
end   

% Case 2
Tm_req_2 = maxTorque-2;  % 78 Nm
Tw_2 = zeros(1,length(omega__wheel));
for i=1:length(omega__wheel)
    Tw_2(i) = MotorModel(Tm_req_2,omega__wheel(i),tau_red,eff_red);
end 

% Case 3
Tm_req_3 = maxTorque-5;  % 75 Nm
Tw_3 = zeros(1,length(omega__wheel));
for i=1:length(omega__wheel)
    Tw_3(i) = MotorModel(Tm_req_3,omega__wheel(i),tau_red,eff_red);
end

% Case 4
Tm_req_4 = maxTorque-10;  % 70 Nm
Tw_4 = zeros(1,length(omega__wheel));
for i=1:length(omega__wheel)
    Tw_4(i) = MotorModel(Tm_req_4,omega__wheel(i),tau_red,eff_red);
end

% Case 5
Tm_req_5 = maxTorque-20;  % 60 Nm
Tw_5 = zeros(1,length(omega__wheel));
for i=1:length(omega__wheel)
    Tw_5(i) = MotorModel(Tm_req_5,omega__wheel(i),tau_red,eff_red);
end

% Case 6
Tm_req_6 = maxTorque-30;  % 50 Nm
Tw_6 = zeros(1,length(omega__wheel));
for i=1:length(omega__wheel)
    Tw_6(i) = MotorModel(Tm_req_6,omega__wheel(i),tau_red,eff_red);
end

% Case 7
Tm_req_7 = maxTorque-50;  % 30 Nm
Tw_7 = zeros(1,length(omega__wheel));
for i=1:length(omega__wheel)
    Tw_7(i) = MotorModel(Tm_req_7,omega__wheel(i),tau_red,eff_red);
end


%% Plot the curves

figure('Name','Wheel Torque','NumberTitle','off'), clf
plot(omega__wheel,Tw_1,'Color',color('blue'),'LineWidth',2);
hold on
plot(omega__wheel,Tw_2,'Color',color('red'),'LineWidth',2);
plot(omega__wheel,Tw_3,'Color',color('green'),'LineWidth',2);
plot(omega__wheel,Tw_4,'Color',color('orange'),'LineWidth',2);
plot(omega__wheel,Tw_5,'Color',color('purple'),'LineWidth',2);
plot(omega__wheel,Tw_6,'Color',color('deepsky_blue'),'LineWidth',2);
plot(omega__wheel,Tw_7,'Color',color('gold'),'LineWidth',2);
grid on
set(gca,'fontsize',16)
xlabel('$\omega_{wheel}$ [rad/s]','FontSize',18)
ylabel('Wheel Torque [Nm]','FontSize',18)
title('Wheel torque - Emrax 208 with FSAE Chimera EVO powertrain','FontSize',18)
legend('req motor torque = 80 Nm','req motor torque = 78 Nm','req motor torque = 75 Nm',...
    'req motor torque = 70 Nm','req motor torque = 60 Nm','req motor torque = 50 Nm',...
    'req motor torque = 30 Nm','location','southwest')
