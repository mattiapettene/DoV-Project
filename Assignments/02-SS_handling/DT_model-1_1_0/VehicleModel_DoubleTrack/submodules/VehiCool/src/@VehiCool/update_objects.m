%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function update_objects( objects, simulation_step )
    % Update the objects of the scenario.
    %
    % Arguments
    % ---------
    %  - objects         -> cell array of objects to update.
    %  - simulation_step -> current simulation step. It must be an integer
    %                       greater than or equal to 1.
    %
    % Usage
    % -----
    %  - VehiCool.update_objects( objects )
    %

    % Parse the inputs
    p = inputParser;
    addRequired( p, 'objects', @iscell );
    addRequired( p, 'simulation_step', @(x) floor(x) == x && x >= 1 );
    parse( p, objects, simulation_step );

    % Update the objects
    for i = 1:length( p.Results.objects )
        % Update the object if the simulation step is a multiple of the
        % deferred update
        if mod( p.Results.simulation_step, ...
                p.Results.objects{i}.deferred_update ) == 0
            % Update the root object
            p.Results.objects{i}.update();
        end

        % Update its children (there is no need to check if the simulation step
        % is a multiple of the deferred update, since we are recursively
        % calling this function)
        if ~isempty( p.Results.objects{i}.children )
            VehiCool.update_objects( p.Results.objects{i}.children, ...
                                     p.Results.simulation_step );
        end
    end

end