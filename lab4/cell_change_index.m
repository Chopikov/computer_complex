function [elem] = cell_change_index(elem)

N = length(elem.points);

for i=1:N
    x(i) = elem.points(i).x;
    y(i) = elem.points(i).y;
end

[min_x, index] = min(x);

result = [];
%result = elem.points(index:N)
for i = index:N
   result =  [result, elem.points(i)];
end

for i = 1:index - 1
   result =  [result, elem.points(i)];
end
elem.points = result;

end