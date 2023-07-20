%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

function out = get_state( obj )
    % Get the current state of the object.
    %
    % Outputs
    % -------
    %  - out -> the next state of the object.
    %
    % Usage
    % -----
    %  - out = obj.get_state()
    %

    % If the state is composed of a single sample, return it
    if size( obj.state, 1 ) == 1
        out = obj.state;
    else
        out = obj.state( obj.idx, : );
    end

end