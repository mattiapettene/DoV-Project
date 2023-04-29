<<<<<<< HEAD
function res = resid_pure_Mz_varFz(P,MZ,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    
    tmp_tyre_data.qBz2 = P(1); 
    tmp_tyre_data.qBz3 = P(2);
    tmp_tyre_data.qDz2 = P(3);
    tmp_tyre_data.qDz7 = P(4);
    tmp_tyre_data.qEz2 = P(5);
    tmp_tyre_data.qEz3 = P(6);
    tmp_tyre_data.qHz2 = P(7);
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       mz0  = MF96_MZ0( 0, ALPHA(i), GAMMA, FZ(i), tmp_tyre_data);
       res = res+(mz0-MZ(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(MZ.^2);

end
=======
function res = resid_pure_Mz_varFz(P,MZ,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    
    tmp_tyre_data.qBz2 = P(1); 
    tmp_tyre_data.qBz3 = P(2);
    tmp_tyre_data.qDz2 = P(3);
    tmp_tyre_data.qDz7 = P(4);
    tmp_tyre_data.qEz2 = P(5);
    tmp_tyre_data.qEz3 = P(6);
    tmp_tyre_data.qHz2 = P(7);
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       mz0  = MF96_MZ0( 0, ALPHA(i), GAMMA, FZ(i), tmp_tyre_data);
       res = res+(mz0-MZ(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(MZ.^2);

end
>>>>>>> ass1_Nats
