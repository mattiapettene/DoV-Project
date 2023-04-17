classdef (StrictDefaults)setControlTargets < matlab.System
    
    %% Public, tunable properties
    properties
        lateralController_type  = 0;
        purePursuit_lookAhead   = 0;
        clothoidBased_lookAhead = 0;
        previewPoint_lookAhead  = 0;
        refRoute_points         = 0;
        road_condition          = 0;
        Lf                      = 0;
    end

    %% Private properties
    properties(Access = private)
        % Initialize a clothoid list for the vehicle route 
        vehRoute = ClothoidList();
        length_vehRoute = 0;
    end

    %% Methods
    methods(Access = protected)
        
    %% setupImpl - used for initialization
    function setupImpl(obj)
        % Build the clothoid list for the vehicle route 
        numOfClothoids_route = size(obj.refRoute_points,1);
        for i = 1:numOfClothoids_route-1
            obj.vehRoute.push_back_G1(obj.refRoute_points(i,1),obj.refRoute_points(i,2),obj.refRoute_points(i,3), obj.refRoute_points(i+1,1),obj.refRoute_points(i+1,2),obj.refRoute_points(i+1,3)); 
        end
        obj.length_vehRoute = obj.vehRoute.length;
    end

    %% stepImpl - called at each time step
    function [targetPoint_latControl,speed_req,road_condition,endOfCircuit] = stepImpl(obj,vehPose,req_speed)
        x_vehCoM = vehPose(1);  % Current vehicle CoM x coord [m]
        y_vehCoM = vehPose(2);  % Current vehicle CoM y coord [m]
        theta_vehCoM = vehPose(3);  % Current vehicle attitude [rad]
        
        road_condition = obj.road_condition;

        switch(obj.lateralController_type)
            case {1,2}
                % Use pure pursuit arc path following or clothoid-based lateral controller.
                % Curvilinear coord of the point belonging to the road middle line
                % that is closest w.r.t. the current vehicle position
                [~,~,s_closest,~,~,~] = obj.vehRoute.closestPoint(x_vehCoM,y_vehCoM);
                if (obj.lateralController_type==1)
                    curvAbscissa_lookAhead = s_closest+obj.purePursuit_lookAhead;
                else
                    curvAbscissa_lookAhead = s_closest+obj.clothoidBased_lookAhead;
                end
                if (curvAbscissa_lookAhead>obj.length_vehRoute)
                    curvAbscissa_lookAhead = obj.length_vehRoute;
                end
                [x_lookAhead,y_lookAhead,theta_lookAhead,curv_lookAhead] = obj.vehRoute.evaluate(curvAbscissa_lookAhead);
                targetPoint_latControl = [x_lookAhead,y_lookAhead,theta_lookAhead,curv_lookAhead];
            case {3,4}
                % Use Stanley lateral controller.  
                % Point of the vehicle in correspondence with the front axle
                x_F = x_vehCoM + obj.Lf*cos(theta_vehCoM);
                y_F = y_vehCoM + obj.Lf*sin(theta_vehCoM);
                % Curvilinear coord of the point belonging to the road middle line
                % that is closest w.r.t. the current vehicle position
                [x_closest,y_closest,s_closest,~,~,~] = obj.vehRoute.closestPoint(x_F,y_F);
                [~,~,theta_closest,curv_closest] = obj.vehRoute.evaluate(s_closest);
                targetPoint_latControl = [x_closest,y_closest,rad2deg(theta_closest),curv_closest];
            case 5
                x_F = x_vehCoM + obj.previewPoint_lookAhead*cos(theta_vehCoM);
                y_F = y_vehCoM + obj.previewPoint_lookAhead*sin(theta_vehCoM);
                [x_closest,y_closest,s_closest,~,~,~] = obj.vehRoute.closestPoint(x_F,y_F);
                [~,~,theta_closest,curv_closest] = obj.vehRoute.evaluate(s_closest);
                ep = sqrt((x_F-x_closest)^2 + (y_F-y_closest)^2);
                targetPoint_latControl = [-ep*sign(atan2((y_F-y_closest),(x_F-x_closest))),0,theta_closest-theta_vehCoM,0];
        end

        
        % desired vehicle speed [m/s]
        speed_req = req_speed/3.6;  
%         if (x_vehCoM>=obj.x_parkLot && y_vehCoM>obj.y_lowBound_parkLot && y_vehCoM<obj.y_uppBound_parkLot)   % (s_closest >= obj.vehRoute.length-10)
%             speed_req = 0.1;  % brake when reaching the parking lot
%             endOfCircuit = 1; % flag to indicate whether the end of the circuit has been reached or not
%         else
            endOfCircuit = 0;
%         end
    end

    function [sz_1,sz_2,sz_3,sz_4] = getOutputSizeImpl(~) 
        sz_1 = [1 4];
        sz_2 = [1];
        sz_3 = [1 4];
        sz_4 = [1];
    end

    function [fz1,fz2,fz3,fz4] = isOutputFixedSizeImpl(~)
        fz1 = true;
        fz2 = true;
        fz3 = true;
        fz4 = true;
    end

    function [dt1,dt2,dt3,dt4] = getOutputDataTypeImpl(~)
        dt1 = 'double';
        dt2 = 'double';
        dt3 = 'double';
        dt4 = 'double';
    end

    function [cp1,cp2,cp3,cp4] = isOutputComplexImpl(~)
        cp1 = false;
        cp2 = false;
        cp3 = false;
        cp4 = false;
    end
        
    end
end
