function r = get_curv_radius(A, B, C)
    a = B - A;
    b = C - B;
    c = C - A;
    cos = sum(c.*b)/(norm(c)*norm(b));
    sin = sqrt(1 - cos^2);
    r = norm(a) / (2 *sin);
end