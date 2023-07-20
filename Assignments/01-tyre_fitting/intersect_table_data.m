function [out_table, out_pos] = intersect_table_data(varargin)
% Intersect tables sequentially

  % Check input tables
  if isrow(varargin{1})
    error('intersectTableData(): Please provide more than 1 table!')
  end

  % Intersect tables sequentially
  [out_table, out_pos] = intersect(varargin{1}, varargin{2}, 'stable');
  if nargin > 2
    for i = 3:nargin
      [out_table, out_pos] = intersect(out_table, varargin{i}, 'stable');
    end
  end
  
  if height(out_table) == 0
    out_table{1,:} = NaN;
  end
  
end

