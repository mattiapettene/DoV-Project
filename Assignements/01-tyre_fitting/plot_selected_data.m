function plot_selected_data(tab)

  color     = '#0072BD';
  linewidth = 1.1;
  linestyle = 'none';
  marker    = '*';

  % Plot selected data
   tiledlayout(3,2);

  % FX ---------
  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\kappa$ (-)');
  ylabel('$F_x$ (N)');
  plot(tab.SL, tab.FX, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\alpha$ (deg)');
  ylabel('$F_x$ (N)');
  plot(tab.SA, tab.FX, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;
 xlim([-0.5 0.5])
  % FY ---------
  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\kappa$ (-)');
  ylabel('$F_y$ (N)');
  plot(tab.SL, tab.FY, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\alpha$ (deg)');
  ylabel('$F_y$ (N)');
  plot(tab.SA, tab.FY, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;
  xlim([-0.5 0.5])
  % MZ ---------
  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\kappa$ (-)');
  ylabel('$M_z$ (Nm)');
  plot(tab.SL, tab.MZ, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\alpha$ (deg)');
  ylabel('$M_z$ (Nm)');
  plot(tab.SA, tab.MZ, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;
 xlim([-0.5 0.5])
  hold off;
end