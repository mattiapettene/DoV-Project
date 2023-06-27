%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function update( obj )
    % Update the object's position and orientation
    %
    % Usage
    % -----
    %  - obj.update()
    %

    % Get the current state
    state = obj.get_state();

    % Create the transform matrix
    T = makehgtform( 'translate', [state(1), state(2), state(3)], ...
                     obj.rot_seq{1}, state(4),     ...
                     obj.rot_seq{2}, state(5),     ...
                     obj.rot_seq{3}, state(6) );

    % Apply the transform
    obj.transform( T );

    % Update the bounding box
    obj.update_bounding_box();

    % Increment the index
    obj.idx = obj.idx + 1;

end