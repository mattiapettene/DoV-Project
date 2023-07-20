function VerticalLoads = LoadTransferModel(Omega,u,u__dot,FAz__r,FAz__f,Mz_tires,params)

    % ----------------------------------------------------------------
    %% Function purpose: compute vertical loads for the 4 wheels 
    % ----------------------------------------------------------------
    VerticalLoads.Fz_rr_ss = 0;
    VerticalLoads.Fz_rl_ss = 0;
    VerticalLoads.Fz_fr_ss = 0;
    VerticalLoads.Fz_fl_ss = 0;

    % Extract the vehicle data
    Ks__r = params.rear_suspension.Ks_r;
    Ks__f = params.front_suspension.Ks_f;
    hrr   = params.rear_suspension.h_rc_r;
    hrf   = params.front_suspension.h_rc_f;
    hG    = params.vehicle.hGs;
    Lr = params.vehicle.Lr;
    Lf = params.vehicle.Lf;
    Wr = params.vehicle.Wr;
    Wf = params.vehicle.Wf;
    m  = params.vehicle.m;
    g  = params.vehicle.g;
    hr = hrf + (hrr-hrf)*Lf/(Lr+Lf);  % [m] height from ground of the projection of the CoM G on the roll axis
    hs = hG - hr;
    
    % Extract the tire self-aligning torques
    AM__rr = Mz_tires.rr;
    AM__rl = Mz_tires.rl;
    AM__fr = Mz_tires.fr;
    AM__fl = Mz_tires.fl;

    % Compute the steady-state vertical loads
    Fz__rr_ss = (((Lf + Lr) * Wr * (Ks__f + Ks__r) * FAz__r + 2 * ((hrr * Lf + hs * (Lf + Lr)) * Ks__r + hrr * Ks__f * Lf) * m * Omega * u + (2 * hrr * AM__fl + 2 * hrr * AM__fr + 2 * hrr * AM__rl + 2 * hrr * AM__rr + Wr * m * (g * Lf + hG * u__dot)) * (Ks__f + Ks__r)) / (Lf + Lr) / Wr / (Ks__f + Ks__r)) / 0.2e1;
    Fz__rl_ss = (((Lf + Lr) * Wr * (Ks__f + Ks__r) * FAz__r - 2 * ((hrr * Lf + hs * (Lf + Lr)) * Ks__r + hrr * Ks__f * Lf) * m * Omega * u + (Ks__f + Ks__r) * (-2 * hrr * AM__fl - 2 * hrr * AM__fr - 2 * hrr * AM__rl - 2 * hrr * AM__rr + Wr * m * (g * Lf + hG * u__dot))) / (Lf + Lr) / Wr / (Ks__f + Ks__r)) / 0.2e1;
    Fz__fr_ss = (((Lf + Lr) * Wf * (Ks__f + Ks__r) * FAz__f + 2 * m * ((hrf * Lr + hs * (Lf + Lr)) * Ks__f + hrf * Ks__r * Lr) * Omega * u + (-2 * hrf * AM__fl - 2 * hrf * AM__fr - 2 * hrf * AM__rl - 2 * hrf * AM__rr + Wf * m * (g * Lr - hG * u__dot)) * (Ks__f + Ks__r)) / (Lf + Lr) / Wf / (Ks__f + Ks__r)) / 0.2e1;
    Fz__fl_ss = (((Lf + Lr) * Wf * (Ks__f + Ks__r) * FAz__f - 2 * m * ((hrf * Lr + hs * (Lf + Lr)) * Ks__f + hrf * Ks__r * Lr) * Omega * u + (Ks__f + Ks__r) * (2 * hrf * AM__fl + 2 * hrf * AM__fr + 2 * hrf * AM__rl + 2 * hrf * AM__rr + Wf * m * (g * Lr - hG * u__dot))) / (Lf + Lr) / Wf / (Ks__f + Ks__r)) / 0.2e1;

    if Fz__rr_ss>=0
        VerticalLoads.Fz_rr_ss = Fz__rr_ss;
    end

    if Fz__rl_ss>=0
        VerticalLoads.Fz_rl_ss = Fz__rl_ss;
    end

    if Fz__fr_ss>=0
        VerticalLoads.Fz_fr_ss = Fz__fr_ss;
    end

    if Fz__fl_ss>=0
        VerticalLoads.Fz_fl_ss = Fz__fl_ss;
    end

end

