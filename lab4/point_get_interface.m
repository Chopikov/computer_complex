function [result] = point_get_interface()
    result = struct('create', @point_create, 'dist', @point_dist);
end
