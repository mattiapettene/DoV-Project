%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function plot( obj, ax )
    % This function sets the camera.
    %
    % Arguments
    % ----------
    %  - ax -> axis handle.
    %
    % Usage
    % -----
    %  - obj.plot( ax )
    %

    % Parse the inputs
    p = inputParser();
    p.addRequired( 'ax', @ishandle );
    p.parse( ax );

    % Set axis handle
    obj.ax = p.Results.ax;

    % Set camera projection
    camproj( obj.ax, obj.projection );

    % Update the camera
    obj.update();

end