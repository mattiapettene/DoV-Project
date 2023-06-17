function Tw = brakeModelRear(reqBrakeTorque,omega,params)

    % ----------------------------------------------------------------
    %% Function purpose: compute braking torques at rear wheels with brake model
    % ----------------------------------------------------------------
    
    max_brake_torque_rear = params.braking.max_brake_torque_rear;
    regularSignScale = params.braking.regularSignScale;
    
    % Check that the braking torques have correctly been specified as
    % negative quantities. Otherwise, set them to zero
    if (reqBrakeTorque > 0)
        reqBrakeTorque = 0;
    end
    
    % Check that the requested braking torque is lower than the one that
    % the hydraulic system can apply. 
    % To make sure that the 
    % vehicle correctly stops at zero forward speed and to avoid that negative 
    % speed values could be reached during braking, braking torque must
    % become zero when the wheel angular speed is negative or zero
    if (abs(reqBrakeTorque) < abs(max_brake_torque_rear))
        Tw = reqBrakeTorque*sin(atan(omega/regularSignScale));
    else
        Tw = -max_brake_torque_rear*sin(atan(omega/regularSignScale));
    end

end

