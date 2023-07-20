function [gamma__rr,gamma__rl,gamma__fr,gamma__fl] = camberModel(params)

    % ----------------------------------------------------------------
    %% Function purpose: compute wheel camber values
    % ----------------------------------------------------------------
    
    rear_static_camber  = params.rear_wheel.static_camber;
    front_static_camber = params.front_wheel.static_camber;
    
    gamma__rr = -rear_static_camber*pi/180;    
    gamma__rl = rear_static_camber*pi/180;   
    gamma__fr = -front_static_camber*pi/180;
    gamma__fl = front_static_camber*pi/180;
    
end

