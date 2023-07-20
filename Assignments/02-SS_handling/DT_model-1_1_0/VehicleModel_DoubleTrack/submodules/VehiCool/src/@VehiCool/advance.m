%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function advance( obj, simulation_step )
    % Advance the scenario by one step.
    %
    % Arguments
    % ---------
    %  - simulation_step -> the current simulation step. It must be an integer
    %                       greater than or equal to 1.
    %
    % Usage
    % -----
    %  - obj.advance()
    %

    % Parse the inputs
    p = inputParser;
    p.addRequired( 'simulation_step', @(x) floor(x) == x && x >= 1 );
    p.parse( simulation_step );

    % Update the objects
    VehiCool.update_objects( obj.objects, p.Results.simulation_step );

    % Update the camera
    obj.camera.update();

end