function [FAxc,FAz__r,FAz__f] = aeroModel(u,vehicle_data)

    % ----------------------------------------------------------------
    %% Function purpose: compute aerodynamic drag and downforce on front 
    %%                   and rear axles
    % ----------------------------------------------------------------

    % Load vehicle aero data
    CAx  = vehicle_data.aerodynamics.CAx;   
    CAzf = vehicle_data.aerodynamics.CAzf;  
    CAzr = vehicle_data.aerodynamics.CAzr;

    FAxc   = CAx*u^2;
    FAz__r = CAzr*u^2;
    FAz__f = CAzf*u^2;

end

