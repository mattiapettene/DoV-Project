function res = resid_pure_Fx(P,FX,KAPPA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
       
    tmp_tyre_data.pCx1 = P(1); 
    tmp_tyre_data.pDx1 = P(2);
    tmp_tyre_data.pEx1 = P(3);
    tmp_tyre_data.pEx4 = P(4);
    tmp_tyre_data.pHx1 = P(5);
    tmp_tyre_data.pKx1 = P(6);
    tmp_tyre_data.pVx1 = P(7);
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(KAPPA)
       fx0  = MF96_FX0(KAPPA(i), 0, GAMMA, FZ, tmp_tyre_data);
       res = res+(fx0-FX(i))^2;
       %res = res+(fx0/FX(i)-1)^2;
    end
    
    % Compute the residuals
    res = res/sum(FX.^2);
    

end

