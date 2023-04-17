function [kappa__rr_corr,kappa__rl_corr,kappa__fr_corr,kappa__fl_corr] = longitSlipModel(u,Omega,kappa__rr,kappa__rl,kappa__fr,kappa__fl,omega__rr,omega__rl,omega__fr,omega__fl,delta__fr,delta__fl,Fz__rr,Fz__rl,Fz__fr,Fz__fl,gamma__rr,gamma__rl,gamma__fr,gamma__fl,pacejkaParam,vehicle_data)

    % ----------------------------------------------------------------
    %% Function purpose: compute low-speed corrections for longitudinal slips
    % ----------------------------------------------------------------
    
    Rf = vehicle_data.front_wheel.Rf;
    Rr = vehicle_data.rear_wheel.Rr;
    Wf = vehicle_data.vehicle.Wf;    
    Wr = vehicle_data.vehicle.Wr; 
    Lf = vehicle_data.vehicle.Lf;        
    Vlow = vehicle_data.tire.Vlow_long;

    % Load Pacejka fitted coeffs
    pCx1 = pacejkaParam.pCx1;
    pDx1 = pacejkaParam.pDx1;
    pDx2 = pacejkaParam.pDx2;
    pDx3 = pacejkaParam.pDx3;
    pKx1 = pacejkaParam.pKx1;
    
    pKx2 = pacejkaParam.pKx2;
    pKx3 = pacejkaParam.pKx3;

    Fz0 = pacejkaParam.Fz0;

    % Pacejka equations
    Fz01 = Fz0;
    % Perform the calculations for each of the 4 wheels
    dfz_rr = Fz__rr / Fz01 - 1;
    dfz_rl = Fz__rl / Fz01 - 1;
    dfz_fr = Fz__fr / Fz01 - 1;
    dfz_fl = Fz__fl / Fz01 - 1;

    Cx = pCx1;
    
    % rear right wheel
    mu__x_rr = (dfz_rr * pDx2 + pDx1) * (-pDx3 * gamma__rr ^ 2 + 1);
    Dx_rr = mu__x_rr * Fz__rr;
    Kxk_rr = Fz__rr * (dfz_rr * pKx2 + pKx1) * exp(-(pKx3 * dfz_rr));
    Bx_rr = Kxk_rr / Cx / Dx_rr;
    CFk_rr = Bx_rr*Cx*Dx_rr;
    
    % rear left wheel
    mu__x_rl = (dfz_rl * pDx2 + pDx1) * (-pDx3 * gamma__rl ^ 2 + 1);
    Dx_rl = mu__x_rl * Fz__rl;
    Kxk_rl = Fz__rl * (dfz_rl * pKx2 + pKx1) * exp(-(pKx3 * dfz_rl));
    Bx_rl = Kxk_rl / Cx / Dx_rl;
    CFk_rl = Bx_rl*Cx*Dx_rl;
    
    % front right wheel
    mu__x_fr = (dfz_fr * pDx2 + pDx1) * (-pDx3 * gamma__fr ^ 2 + 1);
    Dx_fr = mu__x_fr * Fz__fr;
    Kxk_fr = Fz__fr * (dfz_fr * pKx2 + pKx1) * exp(-(pKx3 * dfz_fr));
    Bx_fr = Kxk_fr / Cx / Dx_fr;
    CFk_fr = Bx_fr*Cx*Dx_fr;
    
    % front left wheel
    mu__x_fl = (dfz_fl * pDx2 + pDx1) * (-pDx3 * gamma__fl ^ 2 + 1);
    Dx_fl = mu__x_fl * Fz__fl;
    Kxk_fl = Fz__fl * (dfz_fl * pKx2 + pKx1) * exp(-(pKx3 * dfz_fl));
    Bx_fl = Kxk_fl / Cx / Dx_fl;
    CFk_fl = Bx_fl*Cx*Dx_fl;
    
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

