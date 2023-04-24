function res = resid_pure_Fy_varFz(P,FY,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    
    tmp_tyre_data.pDy2 = P(1); 
    tmp_tyre_data.pEy2 = P(2);
    tmp_tyre_data.pHy2 = P(3);
    tmp_tyre_data.pVy2 = P(4);
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       fy0  = MF96_FY0( 0, ALPHA(i), GAMMA, FZ(i), tmp_tyre_data);
       res = res+(fy0-FY(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FY.^2);

end
