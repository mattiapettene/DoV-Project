%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function plot( obj, ax )
    % Plot the object
    %
    % Arguments
    % ---------
    %  - ax -> the axes handle.
    %
    % Usage
    % -----
    %  - obj.plot( ax )
    %

    % Set the axes handle
    obj.ax = ax;

    % Create the object descriptor
    obj.descriptor = obj.create_descriptor();

    % Create the transform object
    obj.patch_transform = hgtransform( 'Parent', obj.ax );

    % Set the patch object to be a child of the transform object
    set( obj.descriptor, 'Parent', obj.patch_transform );

    % Set the initial transform
    set( obj.patch_transform, 'Matrix', obj.default_transform );

    % Update the object's position
    obj.update();

end