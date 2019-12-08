function [result, count] = get_intersection(elem, H)
POINT = point_get_interface();
ELEMENT = cell_get_interface();
N = length(elem.points);
    
result_left=[];
result_right=[];
N_right = 0;
N_left = 0; 
for i = 1:N
    z = elem.points(i).y;
    r = elem.points(i).x;
    y1 = sqrt(r*r-H*H);
    y2 = -sqrt(r*r-H*H);
    
    if(isreal(y1))
        N_right = N_right+1;
        new_point = POINT.create(y1, z);
        result_right = [result_right, new_point];

    end
    if(isreal(y2))
        N_left = N_left+1;
        new_point = POINT.create(y2, z);
        result_left = [result_left, new_point];
    end
end

result_left = rot90(result_left, 2);
% ни одна из точек не попала на сечение
result = [];
if(N_right == 0)
    count = 0;
end

% не все точки присутствуют на сечении, правый и левый элементы склеиваем
if(0  < N_right && N_right < N)
    count = 1;
    new_elem = ELEMENT.create();
    new_elem.points = [result_right, result_left];
    new_elem.index = elem.index;
    result = [new_elem];
end

% все точки присутствуют на сечении, поличили 2 элемента

if(N_right == N)
    count = 2;
    new_elem1 = ELEMENT.create();
    new_elem1.points = [result_right];
    new_elem1.index = elem.index;
    new_elem2 = ELEMENT.create();
    new_elem2.points = [result_left];
    new_elem2.index = elem.index;
    result = [new_elem1, new_elem2];
end





        