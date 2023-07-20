function res = resid_Fy_varAlpha(P,FY,KAPPA,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    % assigne computed parameter
    tmp_tyre_data.rBy1 = P(1);  
    tmp_tyre_data.rBy2 = P(2); 
    tmp_tyre_data.rBy3 = P(3); 
    tmp_tyre_data.rCy1 = P(4);
    tmp_tyre_data.rHy1 = P(5);
    tmp_tyre_data.rVy1 = P(6);
    tmp_tyre_data.rVy4 = P(7);
    tmp_tyre_data.rVy5 = P(8);
    tmp_tyre_data.rVy6 = P(9);
        
    % Longitudinal Force Equations
    res = 0;
    for i=1:length(KAPPA)
       [fy_fit,~]  = MF96_FY(KAPPA(i), ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data);
       res = res+(fy_fit-FY(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FY.^2);

end

