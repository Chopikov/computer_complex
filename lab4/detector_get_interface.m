function [result] = detector_get_interface()
    result = struct('get_ray', @detector_get_ray, 'get_plane', @detector_get_plane, 'get_col_pos', @detector_get_col_pos);   
end