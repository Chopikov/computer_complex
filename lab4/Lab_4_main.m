%%
clear all
POINT = point_get_interface();
ELEMENT = cell_get_interface();
DETECTOR = detector_get_interface();

%%
[flux, RBDRY, ZBDRY, NBDRY, R, Z, time, rdim, zdim] = gfile_extractor_1t(37000, 000155, 65);

%обход против часовой 
[RBDRY,ZBDRY] = circle_spin_and_reverse(RBDRY,ZBDRY,NBDRY, 32);

[arr, ind_arr] = min(flux);
[flux_min, min_j] = min (arr);
min_i = ind_arr(min_j);
magnet_axis = [R(min_j), Z(min_i)];

figure()
grid on
hold on
plot(RBDRY, ZBDRY, "o");
plot(magnet_axis(1), magnet_axis(2), "*");
[x_limit, y_limit] = get_plot_lim(RBDRY, ZBDRY);
xlim(x_limit) 
ylim(y_limit) 
axis equal
title("separatrix")
legend("separatrix", "magnetic axis")

%%
num_sectors = 8;
num_circle = 6;
cur_cut_y = RBDRY;
cur_cut_z = ZBDRY;
cur_magnet_axis = magnet_axis;

[R_segments_arr, Z_segments_arr, lines_start, lines_end] = get_separatrix_grid(cur_cut_y, cur_cut_z, cur_magnet_axis, num_sectors, num_circle);
elements = create_cells(R_segments_arr, Z_segments_arr, lines_start, lines_end);

figure()
grid on
hold on
axis equal
N = length(elements);
for i = 1:N
    elem = elements(i)
    len = length(elem.points);
    x = zeros(len+1, 1);
    y = zeros(len+1, 1);
    for j=1:len
        x(j) = elem.points(j).x;
        y(j) = elem.points(j).y;
    end

    txt_pos_x = sum(x) / (len);
    txt_pos_y = sum(y) / (len);

    x(len+1) = elem.points(1).x;
    y(len+1) = elem.points(1).y;

    if i < 20 || i==N
        plot(x, y, "b");
        str_num = "";
        text(txt_pos_x, txt_pos_y, strcat(num2str(elem.index),str_num));
    elseif(i == 20)
        plot(x, y, "b");
        text(txt_pos_x, txt_pos_y, "...");
    else
        plot(x, y, "b");
    end
end

%% данные о детекторе
ang = acos((708^2 + 720^2 - 31^2) / (2 * 708 * 720));
spd_start = [0, -0.708];
spd_end = [0.72 * sin(ang), 0.72 * -cos(ang)];
spd_vect = (spd_end - spd_start) / norm(spd_end - spd_start);
spd_xy_step = [2.3375 - 0.88 , 3.81 - 2.3375 + 0.88 ] * 1e-03;
pp = spd_start + spd_vect * ((spd_xy_step(1) + spd_xy_step(2)) * 8 + 0.52 * 1e-03) / 2;
aperture_xy_offset = 0.0395;
aperture_xy = [pp(1) - spd_vect(2) * aperture_xy_offset, pp(2) + spd_vect(1) * aperture_xy_offset];
spd_z_start = (27.52 - 0.49) / 2 * 1e-03;
spd_z_step = -1.72 * 1e-03;
spd_xy = spd_start + spd_vect * (spd_xy_step(2) / 2 + 0.26 * 1e-03);
%%
p = POINT.create(0, 0);
detector = struct('start', p, 'end', p,'step', [], 'direction', p, 'aperture_offset', 0, 'center', p, 'aperture_pos', p, 'horizontal_step', p, 'z_start', 0, 'z_step', 0);   
detector.start = POINT.create(spd_start(1), spd_start(2));
detector.end = POINT.create(spd_end(1), spd_end(2));
detector.step = spd_xy_step;
detector.direction = POINT.create(spd_vect(1), spd_vect(2));
detector.center = POINT.create(pp(1), pp(2));
detector.aperture_offset = aperture_xy_offset;
detector.aperture_pos = POINT.create(aperture_xy(1), aperture_xy(2)); 
detector.z_start = spd_z_start;
detector.z_step = spd_z_step;
detector.horizontal_step = POINT.create((detector.end.x - detector.start.x) / 16 , (detector.end.y - detector.start.y) / 16 ); 

figure()
grid on
hold on
axis equal
phi = -pi:pi/360:pi;
R = max(RBDRY);
r = min(RBDRY);
plot(R * cos(phi), R * sin(phi));
plot(r * cos(phi), r * sin(phi));
plot(detector.aperture_pos.x, detector.aperture_pos.y, "ob")

x = detector.start.x;
y = detector.start.y;
for i=0:16
    plot(x + i * detector.horizontal_step.x, y + i * detector.horizontal_step.y, "or")  
end
plot(detector.center.x, detector.center.y, "*b");
plot(detector.start.x, detector.start.y, "*b") 
plot(detector.end.x, detector.end.y, "*b")

for i=1:16
    A = DETECTOR.get_col_pos(detector, i);
    x = A.x:0.01:0.8;
    B = detector.aperture_pos;
    [k, b] = get_line(A, B);
    plot(x, k*x+b, "r")
    if i == 1 || i == 16
        text_R = R*1.2;
        text_x = (-2*b*k + sqrt(4*b^2*k^2 - 4*(k^2+1)*(b^2 - text_R^2)))/(2*(k^2+1));
        text_y = k*text_x + b;
        text(text_x, text_y, num2str(i));
    end
end

xlim([-0.8, 0.8])
ylim([-0.8, 0.8])

%%
min_x_sep = min(RBDRY);
max_x_sep = max(RBDRY);
cut_ind = 16;
figure()
grid on
hold on
axis equal
H = DETECTOR.get_plane(detector, cut_ind);
title("ind = 16")

h = length(elements);
cut_elements = [];
draw_cut_elements = [];
for i = 1:h
    [elem, count] = get_intersection(elements(i), H);
    draw_cut_elements = [draw_cut_elements, elem];
    if(H < min_x_sep)
        if(count == 2) 
            cut_elements = [cut_elements, elem(2)];
        else
            cut_elements = [cut_elements, elem];
        end
    else
        cut_elements = [cut_elements, elem];
    end
end
N = length(draw_cut_elements);
for i = 1:N
    elem = draw_cut_elements(i)
    len = length(elem.points);
    x = zeros(len+1, 1);
    y = zeros(len+1, 1);
    for j=1:len
        x(j) = elem.points(j).x;
        y(j) = elem.points(j).y;
    end
    x(len+1) = elem.points(1).x;
    y(len+1) = elem.points(1).y;
    plot(x, y, "b");
end
N = length(cut_elements); 
for ray_ind=1:16
    [k, b, det_pos, apper_pos] = DETECTOR.get_ray(detector, cut_ind, ray_ind, true);
    x = det_pos.x:0.01:0.7;
    plot(x, k*x+b, "r")
    plot(det_pos.x, det_pos.y, "*r")

    for t = 1:N
        [hord, intersection] = ELEMENT.get_hord(cut_elements(t), k, b);          
        for g = 1:length(intersection)
            plot(intersection(g).x, intersection(g).y, "om");
        end
    end
    if ray_ind == 1 || ray_ind == 16
        text_x = 0.6;
        text_y = text_x * k + b;
        text(text_x, text_y, num2str(ray_ind));
    end
end
plot(apper_pos.x, apper_pos.y, "ob")
%%
for cut_ind=1:16
    figure()
    grid on
    hold on
    axis equal
    H = DETECTOR.get_plane(detector, cut_ind);
    title1 = strcat("H = ", num2str(H))
    title2 = strcat("(column ", num2str(cut_ind), ")")
    title(title1 + title2);
    element_num = length(elements);
    cut_elements = [];
    for i = 1:element_num
        [elem, count] = get_intersection(elements(i), H);
        if(H < min_x_sep)
            if(count == 2)
                cut_elements = [cut_elements, elem(2)];
            else
                cut_elements = [cut_elements, elem];
            end
        else
            cut_elements = [cut_elements, elem];
        end
    end
    N = length(cut_elements);
    
    for i = 1:N
        elem = cut_elements(i)
        len = length(elem.points);
        x = zeros(len+1, 1);
        y = zeros(len+1, 1);
        for j=1:len
            x(j) = elem.points(j).x;
            y(j) = elem.points(j).y;
        end
        x(len+1) = elem.points(1).x;
        y(len+1) = elem.points(1).y;
        plot(x, y, "b");
    end

    for ray_ind=1:16
        [k, b, det_pos, apper_pos] = DETECTOR.get_ray(detector, cut_ind, ray_ind);
        x = det_pos.x:0.01:0.6;
        plot(x, k*x+b, "r")
        plot(det_pos.x, det_pos.y, "*r")
        for t = 1:N
            [hord, intersection] = ELEMENT.get_hord(cut_elements(t), k, b);          
            for g = 1:length(intersection)
                plot(intersection(g).x, intersection(g).y, "ok");
            end
        end
    end    
    plot(apper_pos.x, apper_pos.y, "ob")
end
%%
%график плоскости детектора
figure()
grid on
hold on
spd_sizes = [0.88, 1.23] * 1e-3;

y_square_up = [0, spd_sizes(1), spd_sizes(1), 0, 0];
y_square_down = [0, 0, spd_sizes(1), spd_sizes(1), 0];
z_square_up = [0, 0, spd_sizes(2), spd_sizes(2), 0];
z_square_down = [0, -spd_sizes(2), -spd_sizes(2), 0, 0];
color = "k";
point_color = ".k";
for j = 0:7
    for i = 0:7
        x = (spd_sizes(1) + j * sum(spd_xy_step) + y_square_up);
        y = (-spd_z_step - spd_sizes(2)) / 2 + i * (-spd_z_step) + z_square_up;
        plot(x, y, color);
        x = x(1:4);
        x = sum(x)/4;
        y = y(1:4);
        y = sum(y)/4;
        plot(x, y, point_color);
        x = spd_sizes(1) + spd_xy_step(1) + j * sum(spd_xy_step) + y_square_up;
        y = (-spd_z_step - spd_sizes(2)) / 2 + i * (-spd_z_step) + z_square_up;
        plot(x, y, color);
        x = x(1:4);
        x = sum(x)/4;
        y = y(1:4);
        y = sum(y)/4;
        plot(x, y, point_color);
        x = spd_sizes(1) + j * sum(spd_xy_step) + y_square_down;
        y = (spd_z_step + spd_sizes(2)) / 2 + i * spd_z_step + z_square_down;
        plot(x, y, color);
        x = x(1:4);
        x = sum(x)/4;
        y = y(1:4);
        y = sum(y)/4;
        plot(x, y, point_color);
        x = spd_sizes(1) + spd_xy_step(1) + j * sum(spd_xy_step) + y_square_down;
        y = (spd_z_step + spd_sizes(2)) / 2 + i * spd_z_step + z_square_down;
        plot(x, y, color);
        x = x(1:4);
        x = sum(x)/4;
        y = y(1:4);
        y = sum(y)/4;
        plot(x, y, point_color);
    end
end
%%
element_num = length(elements);
hord_matrix = zeros(256, element_num);
cur_row = 1;
for cut_ind=1:16
    H = DETECTOR.get_plane(detector, cut_ind);
    cut_elements = [];
    for i = 1:element_num
        [elem, count] = get_intersection(elements(i), H);
        if(H < min_x_sep)
            if(count == 2)
                cut_elements = [cut_elements, elem(2)];
            else
                cut_elements = [cut_elements, elem];
            end
        else
            cut_elements = [cut_elements, elem];
        end
    end
    N = length(cut_elements);
    for ray_ind=1:16
        [k, b, det_pos, apper_pos] = DETECTOR.get_ray(detector, 16, ray_ind);
        for t = 1:N
            [hord, intersection] = ELEMENT.get_hord(cut_elements(t), k, b); 
            hord_matrix(cur_row, cut_elements(t).index) = hord_matrix(cur_row, cut_elements(t).index) + hord;
        end
        cur_row = cur_row + 1;
    end
end
