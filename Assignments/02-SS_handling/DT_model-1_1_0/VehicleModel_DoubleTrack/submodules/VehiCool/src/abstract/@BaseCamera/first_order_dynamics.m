%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function out = first_order_dynamics( obj, state, next_state )
    % This function computes the next of a system following the
    % dynamics of a first order system.
    %
    % Arguments
    % ----------
    %  - state      -> current state of the system.
    %  - next_state -> theoretical next state of the system.
    %
    % Outputs
    % -------
    %  - out -> next of the system following the dynamics of a first
    %           order system.
    %
    % Usage
    % -----
    %  - out = obj.first_order_dynamics( state, next_state )
    %

    % Parse the inputs
    p = inputParser;
    addRequired( p, 'state', @isnumeric );
    addRequired( p, 'next_state', @isnumeric );
    parse( p, state, next_state );

    % Compute the next state
    out = p.Results.state + obj.delay * ( p.Results.next_state - ...
                                          p.Results.state );

end