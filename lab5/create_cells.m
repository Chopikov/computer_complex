function [result] = create_cells(R_segments_arr, Z_segments_arr, lines_start, lines_end)
    ELEMENT = cell_get_interface();
    POINT = point_get_interface();
    
    elems_result = [];
    N = length(lines_start);
    % добавим начала линий к спискам точек
    tmp_R = zeros(1, N);
    tmp_Z = zeros(1, N);
    for i = 1:N
       tmp_R(i) = lines_start{i}(1);
       tmp_Z(i) = lines_start{i}(2);
    end
    
    R_vector = {tmp_R R_segments_arr{:}};
    Z_vector = {tmp_Z Z_segments_arr{:}};
     
    circle_N = length(R_vector);
    index = 0;
    for j = 1:circle_N-1
        cur_z = Z_vector{j};
        cur_r = R_vector{j};
        next_cur_z = Z_vector{j + 1};
        next_cur_r = R_vector{j + 1};
        
        for i = 1:N - 1
            new_elem = ELEMENT.create();
            p1 = POINT.create(cur_r(i), cur_z(i));
            p2 = POINT.create(cur_r(i+1), cur_z(i+1));
            p3 = POINT.create(next_cur_r(i + 1), next_cur_z(i + 1));
            p4 = POINT.create(next_cur_r(i), next_cur_z(i));
            new_elem.points = [p1, p2, p3, p4]; 
            index = index + 1;
            new_elem.index = index;
            elems_result = [elems_result, new_elem];
        end
        
         new_elem = ELEMENT.create();
         p1 = POINT.create(cur_r(N), cur_z(N));
         p2 = POINT.create(cur_r(1), cur_z(1));
         p3 = POINT.create(next_cur_r(1), next_cur_z(1));
         p4 = POINT.create(next_cur_r(N), next_cur_z(N));
         new_elem.points = [p1, p2, p3, p4];
         index = index + 1;
         new_elem.index = index;
         elems_result = [elems_result, new_elem];
    end
    
    cur_z = Z_vector{circle_N};
    cur_r = R_vector{circle_N};

    for i = 1:N - 1
        new_elem = ELEMENT.create();
        p1 = POINT.create(cur_r(i), cur_z(i));
        p2 = POINT.create(cur_r(i+1), cur_z(i+1));
        p3 = POINT.create(lines_end{i}(1), lines_end{i}(2));
        new_elem.points = [p1, p2, p3];   
        index = index + 1;
        new_elem.index = index;
        elems_result = [elems_result, new_elem];
    end
        
    new_elem = ELEMENT.create();
    p1 = POINT.create(cur_r(N), cur_z(N));
    p2 = POINT.create(cur_r(1), cur_z(1));
    p3 = POINT.create(lines_end{1}(1), lines_end{1}(2));
    new_elem.points = [p1, p2, p3];
    index = index + 1;
    new_elem.index = index;
    elems_result = [elems_result, new_elem];
    
    len = length(elems_result);
    for i =1:len
        elems_result(i) = ELEMENT.reindexing(elems_result(i));
        
    end
    result = elems_result;
    
end