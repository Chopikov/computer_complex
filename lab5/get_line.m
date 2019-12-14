function [k, b] = get_line(A, B)
k = (B.y - A.y) / (B.x - A.x);
b = -(B.y - A.y) / (B.x - A.x) * A.x + A.y;
end
