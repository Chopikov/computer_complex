function [hord_dist, intersection] = cell_get_hord(elem, k, b)
POINT = point_get_interface();
ELEMENT = cell_get_interface();
N = length(elem.points);
    
intersection=[];
tmp = [elem.points, elem.points(1)];
for i = 1:N
    A = tmp(i);
    B = tmp(i+1);
       
    if((k * A.x + b - A.y) * (k * B.x + b - B.y) > 0)
        %пересечения нет
        hord_dist = 0;
        continue;
    end
    
    if( (k * A.x + b - A.y) * (k * B.x + b - B.y) == 0)
        if(k * A.x + b - A.y == 0)
            %прошли прямо через точку А
            intersection = [intersection, A];
        end
        %вторую точку проверим как начало следующей стороны
%         if(k * B.x + b - B.y == 0)
%             %прошли прямо через точку B
%             intersection = [intersection, B];
%         end
        continue
    end
    koef = abs((k * A.x + b - A.y)/(k * B.x + b - B.y));
    P = POINT.create((A.x + koef * B.x) / (1 + koef), (A.y + koef * B.y) / (1 + koef));
    intersection = [intersection, P];
end

if(length(intersection) <= 1)
    %пересечений нет (одно не считается)
    hord_dist = 0;
    %intersection = [];
    return 
end
    

%пересечения есть (притом 2 раза)
%FIXME возможен ли случай, когда точек 3 (все три вряд)? 
%тогда учитывать максимальное расстояние среди возможных???
hord_dist =  POINT.dist(intersection(1), intersection(2));
if(length(intersection) > 2)
    ELEMENT.draw(elem, "b", -1, true)
    disp(length(intersection));
end
end