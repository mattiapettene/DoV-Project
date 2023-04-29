function res = resid_pure_Fx_varFz(P,FX,KAPPA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    
       
    tmp_tyre_data.pDx2 = P(1); 
    tmp_tyre_data.pEx2 = P(2);
    tmp_tyre_data.pEx3 = P(3);
    tmp_tyre_data.pHx2 = P(4);
    tmp_tyre_data.pKx2 = P(5);
    tmp_tyre_data.pKx3 = P(6);
    tmp_tyre_data.pVx2 = P(7);
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(KAPPA)
       fx0  = MF96_FX0(KAPPA(i), 0, GAMMA, FZ(i), tmp_tyre_data);
       res = res+(fx0-FX(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FX.^2);

end

