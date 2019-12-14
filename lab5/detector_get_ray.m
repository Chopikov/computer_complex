function [k, b, A, B] = detector_get_ray(det, col, row, is_need_plot)
DETECTOR = detector_get_interface();
POINT = point_get_interface();
H = DETECTOR.get_plane(det, col);

pixel_pos = DETECTOR.get_col_pos(det, col);

detector_dist = POINT.dist(pixel_pos);
aperture_dist = POINT.dist(det.aperture_pos);

detector_x = sqrt(detector_dist^2 - H^2) * sign(pixel_pos.y);
aperture_x = sqrt(aperture_dist^2 - H^2) * sign(det.aperture_pos.y);

A = POINT.create(detector_x, det.z_start + det.z_step * (row-1));

B = POINT.create(aperture_x, 0);

if(nargin == 4)
    if(is_need_plot)
        plot(A.x, A.y, "or") 
        plot(B.x, B.y, "or") 
    end
end

k = (B.y - A.y) / (B.x - A.x);
b = -(B.y - A.y) / (B.x - A.x) * A.x + A.y;
end