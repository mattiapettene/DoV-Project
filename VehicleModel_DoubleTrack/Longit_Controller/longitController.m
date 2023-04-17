function longCTRparam = longitController()

    % ---------------------------------------------------
    %% Load the parameters of the longitudinal controller
    % ---------------------------------------------------
    
    PID_Gain_Scheduling   = readtable('PID_Gain_Scheduling.txt');
    longCTRparam.speed_vals_PID_sched = PID_Gain_Scheduling.speed;
    longCTRparam.Kp_vals  = PID_Gain_Scheduling.Kp;
    longCTRparam.Ki_vals  = PID_Gain_Scheduling.Ki;
    longCTRparam.Kd_vals  = PID_Gain_Scheduling.Kd;
    longCTRparam.N_vals   = PID_Gain_Scheduling.N;
    longCTRparam.Kaw_vals = PID_Gain_Scheduling.Kaw;
    
end

