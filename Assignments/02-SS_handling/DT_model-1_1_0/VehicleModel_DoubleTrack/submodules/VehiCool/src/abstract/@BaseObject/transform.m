%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function transform( obj, T )
    % Transform the object
    %
    % Arguments
    % ---------
    %  - T -> the transformation matrix.
    %
    % Usage
    % -----
    %  - obj.transform( T )
    %

    % Check if the object has a parent
    if ~isempty( obj.parent )
        % Remove the parent's default transform
        R = obj.parent.get_transform() / obj.parent.default_transform;

        % Create the new transformation matrix
        new_transform = R * obj.default_transform * T;

    else

        % Create the new transformation matrix
        new_transform = T * obj.default_transform;

    end

    % Apply the transformation to the object
    set( obj.patch_transform, 'Matrix', new_transform );

end