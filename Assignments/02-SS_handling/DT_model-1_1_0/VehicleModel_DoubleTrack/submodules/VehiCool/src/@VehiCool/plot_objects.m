%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function plot_objects( objects, ax )
    % Plot the objects of the scenario.
    %
    % Arguments
    % ---------
    %  - objects -> cell array of objects.
    %  - ax      -> axes handle.
    %
    % Usage
    % -----
    %  - VehiCool.plot_objects( objects, ax )
    %

    % Parse the inputs
    p = inputParser;
    addRequired( p, 'objects', @iscell );
    addRequired( p, 'ax', @ishandle );
    parse( p, objects, ax );

    % Unravel the tree of objects
    for i = 1:length( p.Results.objects )
        % Plot the root object
        p.Results.objects{i}.plot( p.Results.ax );

        % Plot its children
        if ~isempty( p.Results.objects{i}.children )
            VehiCool.plot_objects( p.Results.objects{i}.children, ...
                                   p.Results.ax );
        end
    end

end