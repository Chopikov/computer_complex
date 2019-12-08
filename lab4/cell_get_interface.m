function [result] = cell_get_interface()
    result = struct('create', @cell_create, 'reindexing', @cell_change_index, 'get_hord', @cell_get_hord);
end