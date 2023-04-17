classdef (StrictDefaults)clothoidBased_lateralController < matlab.System
    
    %% Public, tunable properties
    properties
        lookAhead = 0;
        Kus       = 0;
        L         = 0;
    end

    %% Private properties
    properties(Access = private)
        % Initialize a clothoid for the vehicle path 
        clothoid = ClothoidCurve();
    end

    %% Methods
    methods(Access = protected)
        
    %% setupImpl - used for initialization
    function setupImpl(obj)
        
    end

    %% stepImpl - called at each time step
    function [delta_req] = stepImpl(obj,targetPoint,vehPose,vehSpeed,endOfCircuit)
        % Extract vehicle center of mass position
        x_vehCoM = vehPose(1);      % Current vehicle CoM x coord [m]
        y_vehCoM = vehPose(2);      % Current vehicle CoM y coord [m]
        theta_vehCoM = vehPose(3);  % Current vehicle attitude [rad]
        
        % Extract the parameters of the target point
        x_target = targetPoint(1);      % Target x coord [m]
        y_target = targetPoint(2);      % Target y coord [m]
        theta_target = targetPoint(3);  % Target vehicle attitude [rad]
        
        % Build the clothoid from the current vehicle pose to the target point
        obj.clothoid.build_G1(x_vehCoM,y_vehCoM,theta_vehCoM, x_target,y_target,theta_target);
        initCurv = obj.clothoid.kappaBegin();  % initial clothoid curvature
        
        % Compute the desired steering angle
        if (~endOfCircuit)
            delta_req = initCurv*(obj.L + obj.Kus*vehSpeed^2);
        else
            delta_req = 0;
        end
    end

    function [sz_1] = getOutputSizeImpl(~) 
        sz_1 = [1];
    end

    function [fz1] = isOutputFixedSizeImpl(~)
        fz1 = true;
    end

    function [dt1] = getOutputDataTypeImpl(~)
        dt1 = 'double';
    end

    function [cp1] = isOutputComplexImpl(~)
        cp1 = false;
    end
        
    end
end
