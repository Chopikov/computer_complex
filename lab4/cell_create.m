function [result] = cell_create()
POINT = point_get_interface();
p = POINT.create(0, 0);
result = struct('points', [p], 'index', -1);    
end