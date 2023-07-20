%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function update( obj )
    % This function updates the camera state and target.
    %
    % Usage
    % -----
    %  - obj.update()
    %

    % Check if the camera state has been initialized
    if isempty( obj.state )

        % Initialize the camera state
        obj.state = obj.create_state();

    else

        % Compute the camera next state
        obj.state = obj.first_order_dynamics( obj.state, obj.create_state() );

    end

    % Check if the target state has been initialized
    if isempty( obj.target_state )

        % Initialize the target state
        target_state = obj.get_target_state();
        obj.target_state = target_state(1, 1:3);

    else

        % Compute the target next state
        target_current_state = obj.target_state(1, 1:3);
        target_next_state    = obj.get_target_state();
        obj.target_state     = obj.first_order_dynamics( ...
                                 target_current_state, ...
                                 target_next_state(1, 1:3) );

    end

    % Set camera state
    campos( obj.ax, obj.state );

    % Set camera target
    camtarget( obj.ax, obj.target_state );

    % Set the camera view angle
    camva( obj.ax, obj.view_angle );

end