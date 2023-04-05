%% Plot all the "raw" data

color     = '#0072BD';
linewidth = 1.1;
linestyle = '-';
marker    = 'none';

% -- changing parameters --------------------------------
nsp = 610;

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.KAPPA, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(-)');
title('Slip Ratio');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.ALPHA, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(rad)');
title('Slip Angle');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.RHO, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(m)');
title('Tire Deflection');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.GAMMA, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(rad)');
title('Camber Angle');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.P, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(Pa)');
title('Inflation Pressure');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.VX, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(m/s)');
title('Road Speed');

hold off;

% -- FX, FY, FZ, MX, MY, MZ --------------------------------
figure;
hold on;
nsp = 610;

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.FX, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(N)');
title('Longitudinal force');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.FY, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(N)');
title('Lateral force');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.FZ, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(N)');
title('Vertical load');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.MX, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker);
%xlabel('Sample (-)');
ylabel('(Nm)');
title('Overturning moment');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.MY, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker)
xlabel('Sample (-)');
ylabel('(Nm)');
title('Rolling Resistance');

nsp = nsp+1;
subplot(nsp);
hold on;
grid on;
grid minor;
plot(table.MZ, 'Color', color, 'LineStyle', linestyle, ...
     'LineWidth', linewidth, 'Marker', marker)
xlabel('Sample (-)');
ylabel('(Nm)');
title('Self-aliging torque');

hold off;

