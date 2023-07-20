function [kappa__rr_ss,kappa__rl_ss,kappa__fr_ss,kappa__fl_ss] = steadyStateLongSlips(omega__rr,omega__rl,omega__fr,omega__fl,Omega,u,delta__fr,delta__fl,params)

    % ----------------------------------------------------------------
    %% Function purpose: compute steady-state longitudinal slips
    % ----------------------------------------------------------------
    
    Rr     = params.rear_wheel.Rr;  
    Rf     = params.front_wheel.Rf; 
    Wr     = params.vehicle.Wr;  
    Wf     = params.vehicle.Wf;     
    Lf     = params.vehicle.Lf; 
    Vlow   = params.tire.Vlow_long;
    
    % kappa_rr
    V_cx_rr = Omega*Wr/2 + u;
    V_sx_rr = V_cx_rr - Rr*omega__rr;
    if V_cx_rr > Vlow
        kappa__rr_ss = -V_sx_rr/V_cx_rr;
    else
        % Low-speed case
        kappa__rr_ss = -2*V_sx_rr/( Vlow + (V_cx_rr)^2/Vlow );
    end
    
    % kappa_rl
    V_cx_rl = -Omega*Wr/2 + u;
    V_sx_rl = V_cx_rl - Rr*omega__rl;
    if V_cx_rl > Vlow
        kappa__rl_ss = -V_sx_rl/V_cx_rl;
    else
        % Low-speed case
        kappa__rl_ss = -2*V_sx_rl/( Vlow + (V_cx_rl)^2/Vlow );
    end
    
    % kappa_fr
    V_cx_fr = 1/2*Omega*Wf + u + Lf*Omega*delta__fr;
    V_sx_fr = V_cx_fr - Rf*omega__fr;
    if V_cx_fr > Vlow
        kappa__fr_ss = -V_sx_fr/V_cx_fr;
    else
        % Low-speed case
        kappa__fr_ss = -2*V_sx_fr/( Vlow + (V_cx_fr)^2/Vlow );
    end
    
    % kappa_fl
    V_cx_fl = -1/2*Omega*Wf + u + Lf*Omega*delta__fl;
    V_sx_fl = V_cx_fl - Rf*omega__fl;
    if V_cx_fl > Vlow
        kappa__fl_ss = -V_sx_fl/V_cx_fl;
    else
        % Low-speed case
        kappa__fl_ss = -2*V_sx_fl/( Vlow + (V_cx_fl)^2/Vlow );
    end
    
end

