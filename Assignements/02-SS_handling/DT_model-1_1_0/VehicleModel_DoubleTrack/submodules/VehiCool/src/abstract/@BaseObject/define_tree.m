%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function define_tree( obj, parent )
    % Define the object tree
    %
    % Arguments
    % ---------
    %  - parent -> the parent object.
    %
    % Usage
    % -----
    %  - obj.define_tree( parent )
    %

    % Set the parent
    obj.parent = parent;

    % Add the object to the parent's list of children
    if ~isempty( obj.parent )
        obj.parent.children{end + 1} = obj;
    end

end