function Tw = MotorModel(Tm_request,omega__wheel,params)

    % ----------------------------------------------------------------
    %% Function purpose: compute driving torque at the wheel with motor model
    % ----------------------------------------------------------------

    % Tw is the output wheel torque calculated on the basis of the function
    % input parameters

    tau_red = params.transmission.tau_red;
    eff_red = params.transmission.eff_red;
    
    % Tm_request is the requested motor torque. It can never exceed 80 Nm, 
    % which is the max torque that the motor can provide, at any rotational
    % speed
    max_mot_torque = params.motor.maxTorque;
    
    % Check that the requested motor torque is >=0, otherwise change sign
    if (Tm_request<0)
        Tm_request = -Tm_request;
    end    
    
    % Compute the motor speed corresponding to the current wheel speed
    omega__mot = omega__wheel*tau_red;
    
    % Compute the angular speed of the motor (in rad/s) for which the
    % requested torque will have to be decreased in order to comply with
    % the 80 kW power limit 
    maxPower_split = params.accumulator.maxPower/2;  % max power that can be supplied to a motor
    if (Tm_request<=max_mot_torque && Tm_request>1e-4)
        omega_limit = (maxPower_split*10^3) / Tm_request;
    else
        omega_limit = (maxPower_split*10^3) / max_mot_torque;
    end
    
    maxRotSpeed = params.motor.maxRotSpeed; % max motor speed (limited by DC voltage of battery pack)
    omega_mot_max = maxRotSpeed*pi/30;
    omega_wheel_max = omega_mot_max / tau_red;
    
    % when the motor reaches a rotational speed very close to the max speed,
    % the motor torque starts to decrease a lot. This is used in the
    % simulation in order to achieve the correct value of max vehicle speed
    % by decreasing the motor torque in a linear way, as the motor speed
    % increases beyond the 'speedForTorqueCut'
    speedForTorqueCut = params.motor.speedForTorqueCut;
    omega_motTorqueCut = speedForTorqueCut*pi/30;
    
    % Torque that the motor can deliver when it is at max speed
    torque_mot_at_max_speed = (maxPower_split*10^3) / omega_motTorqueCut;
    
    % Slope of the linear torque decreasing behavior beyond omega_motTorqueCut
    slopeTorqueCut = torque_mot_at_max_speed / (omega_motTorqueCut-omega_mot_max);
    
    % Compute the torque Tw delivered to the wheel 
    if omega__mot <= omega_limit && Tm_request <= max_mot_torque
        if (omega__mot > omega_motTorqueCut + (Tm_request-torque_mot_at_max_speed)/slopeTorqueCut)
            Tw = (torque_mot_at_max_speed + slopeTorqueCut*(omega__mot-omega_motTorqueCut))*tau_red*eff_red;
        else
            Tw = Tm_request*tau_red*eff_red;
        end    
    elseif omega__mot > omega_limit && omega__mot <= omega_motTorqueCut 
        Tw = (maxPower_split*10^3*tau_red*eff_red) / omega__mot;
    elseif (omega__mot > omega_motTorqueCut && omega__mot < omega_mot_max)
        Tw = (torque_mot_at_max_speed + slopeTorqueCut*(omega__mot-omega_motTorqueCut))*tau_red*eff_red;
    elseif omega__mot <= omega_limit && Tm_request > max_mot_torque
        Tw = max_mot_torque*tau_red*eff_red;
    else 
        Tw = torque_mot_at_max_speed*tau_red*eff_red;
    end
  
end

