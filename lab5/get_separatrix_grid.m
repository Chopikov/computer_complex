function [R_segments_arr, Z_segments_arr, lines_start, lines_end] = get_separatrix_grid(cur_cut_y, cur_cut_z, cur_magnet_axis, N, L)
    cur_rad = get_curv_radius_arr( cur_cut_y,  cur_cut_z, length(cur_cut_y));
    cur_rad_smoothed = medfilt1(cur_rad);

    init_sector = 1:length(cur_cut_y);
    
    sector = init_sector;
    r = cur_rad_smoothed;
    separator_ind = fix(length(sector) / 2);  

    first_sector = sector(1:separator_ind);
    second_sector = sector(separator_ind + 1:length(sector));
    [max_val , first_max_ind] = max(r(first_sector));
    [max_val , second_max_ind] = max(r(second_sector));
   
    second_max_ind = second_max_ind + length(first_sector);
    sep_1 = first_max_ind;
    sep_2 = length(first_sector);
    sep_3 = second_max_ind;
    sep_4 = length(sector);
    
    s1 = sector(1 : sep_1);    
    s2 = sector(sep_1 + 1 : sep_2);
    s3 = sector(sep_2 + 1 : sep_3);
    s4 = sector(sep_3 + 1: sep_4);
    
    
    sectors = {s1, s2, s3, s4};
    result = {};
    for i = 1:length(sectors)
        tmp_sector = sectors{i};
        
        res = {};
        cur_step = ceil(length(tmp_sector) / N);
        cur_left = 1;
        cur_right = cur_step;
        reduse = ceil(length(tmp_sector));
        for j = 1:N
            if(reduse == 0)
                break
            end
            tmp = tmp_sector(cur_left:cur_right);
            res = {res{:} tmp};
            reduse = reduse - length(tmp);
            cur_left = cur_right + 1;
            cur_step = ceil(reduse / (N - j));
            cur_right = cur_right + cur_step;  
        end 
        
        result = {result{:} res{:}};
    end
    sectors = result;

    segments = [];
    middle = cur_magnet_axis;
    R_segments_arr = {};
    Z_segments_arr = {}; 
    R_segments = [];
    Z_segments = [];


    for i = 1:length(sectors)
        tmp_sector = sectors{i};
        cur_ind = tmp_sector(1);
        cur_point = [cur_cut_y(cur_ind),cur_cut_z(cur_ind)];
        split_points = {};
        step = (middle - cur_point) / L;
        for k=1:L-1
            split_points = {split_points{:} cur_point + step * (k)};        
        end
        for k = 1:length(split_points)
            cur_split_point = split_points{k};
            R_segments_arr{k}(i) = cur_split_point(1);
            Z_segments_arr{k}(i) = cur_split_point(2);
        end
    end
    
    lines_start = {};
    lines_end = {};
    for i = 1:length(sectors)
        tmp_sector = sectors{i};
        cur_ind = tmp_sector(1);
        cur_point = [cur_cut_y(cur_ind),cur_cut_z(cur_ind)];
        lines_start{i} = cur_point;
        lines_end{i} = middle;
    end
    
end