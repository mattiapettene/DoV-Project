function [kappa__rr_corr,kappa__rl_corr,kappa__fr_corr,kappa__fl_corr] = long_slip_low_speed(u,Omega,kappa__rr,kappa__rl,kappa__fr,kappa__fl,omega__rr,omega__rl,omega__fr,omega__fl,delta__fr,delta__fl,Fz__rr,Fz__rl,Fz__fr,Fz__fl,gamma__rr,gamma__rl,gamma__fr,gamma__fl, vehicle_data)

    % ----------------------------------------------------------------
    %% Function purpose: compute low-speed corrections for longitudinal slips
    % ----------------------------------------------------------------
    
    Rf = vehicle_data.tyre_data_f.R0;
    Rr = vehicle_data.tyre_data_r.R0;
    Wf = vehicle_data.vehicle.Wf;    
    Wr = vehicle_data.vehicle.Wr; 
    Lf = vehicle_data.vehicle.Lf;        
    Vlow = vehicle_data.tire.Vlow_long;
    
    % rear right wheel
    [~, ~, ~, ~, ~, ~, CFk_rr] = MF96_FX0_coeffs(kappa__rr, 0, gamma__rr, Fz__rr, vehicle_data.tyre_data_r);

    
    % rear left wheel
    [~, ~, ~, ~, ~, ~, CFk_rl] = MF96_FX0_coeffs(kappa__rl, 0, gamma__rl, Fz__rl, vehicle_data.tyre_data_r);

    
    % front right wheel
    [~, ~, ~, ~, ~, ~, CFk_fr] = MF96_FX0_coeffs(kappa__fr, 0, gamma__fr, Fz__fr, vehicle_data.tyre_data_f);

    
    % front left wheel
    [~, ~, ~, ~, ~, ~, CFk_fl] = MF96_FX0_coeffs(kappa__fl, 0, gamma__fl, Fz__fl, vehicle_data.tyre_data_f);

    
    kVlow0 = 770;
    
    if abs(u) <= Vlow
        kVlow = 1/2*kVlow0*(1+cos(pi*abs(u)/Vlow));
    else
        kVlow = 0;
    end
   
    Vcx_rr = Omega*Wr/2 + u;
    Vcx_rl = -Omega*Wr/2 + u;
    Vcx_fr = 1/2*Omega*Wf + u + Lf*Omega*delta__fr;
    Vcx_fl = -1/2*Omega*Wf + u + Lf*Omega*delta__fl;
    
    Vsx_rr = Vcx_rr - Rr*omega__rr;
    Vsx_rl = Vcx_rl - Rr*omega__rl;
    Vsx_fr = Vcx_fr - Rf*omega__fr;
    Vsx_fl = Vcx_fl - Rf*omega__fl;

    % New longitudinal slips, with corrections for low speed
    kappa__rr_corr = kappa__rr - kVlow*(Vsx_rr)/CFk_rr;
    kappa__rl_corr = kappa__rl - kVlow*(Vsx_rl)/CFk_rl;
    kappa__fr_corr = kappa__fr - kVlow*(Vsx_fr)/CFk_fr;
    kappa__fl_corr = kappa__fl - kVlow*(Vsx_fl)/CFk_fl;

end

