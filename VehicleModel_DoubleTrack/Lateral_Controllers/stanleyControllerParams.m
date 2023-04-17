function stanleyParams = stanleyControllerParams()

    % ----------------------------------------------------------------
    %% Function purpose: define the Stanley lateral controller parameters
    % ----------------------------------------------------------------
    
    % ---------------
    %% Parameters used by the controller with the KINEMATIC vehicle model
    % ---------------
    
    % max steering angle [deg] at the front wheel of a single track vehicle
    % model. Note that, for control purposes, this angle might be set to a lower
    % value w.r.t. the real max steering angle of the vehicle
    stanleyParams.maxSteerAngle = 6; 
    
    % position gain for the forward motion
    stanleyParams.positGainFW = 0.5; 
    
    % position gain for the reverse motion --> not used
    stanleyParams.positGainREV = 2.5; 
    
        
    % ---------------
    %% Parameters used by the controller with the DYNAMIC vehicle model
    % ---------------
    
    % yaw rate feedback gain
    stanleyParams.yawRateGain = 0.2;
    
    % steering angle feedback gain
    stanleyParams.steerGain = 0.2;
    
    % tire cornering stiffness (slope at the origin of the Fy vs. alpha curve)
    stanleyParams.cornerStiff = 26285;  % [N/rad]

end

