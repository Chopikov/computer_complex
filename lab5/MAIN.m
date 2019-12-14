%%
clear all
POINT = point_get_interface();
ELEMENT = cell_get_interface();
DETECTOR = detector_get_interface();

input_file_index = 37000;
input_file_name = strcat("data\", num2str(37000), "_SPD16x16.mat");
input_time_period = 000155;

%%
[flux, RBDRY, ZBDRY, NBDRY, R, Z, time, rdim, zdim] = gfile_extractor_1t(input_file_index, input_time_period, 65);

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

%spy(hord_matrix)
%%
input_data = load(input_file_name);
sign_bb = input_data.sign_bb(:, :, :);
cnt_meas = size(sign_bb, 3);
tp = cell2mat(input_data.Data(1, 2)) * 1e-3;% - шаг по времени
tz = cell2mat(input_data.Data(2, 2));
t_start = tz;
t_end = t_start + (cnt_meas - 1) * tp;
t_i = t_start:tp:t_end;

%%
% выбираем "окно" по котору мычисляем границы b
b_time_window = 2;

%извлекаем значения по нужным периода времени 
b_data = [];
for start_efit_time1 = input_time_period - b_time_window:input_time_period + b_time_window
    ind = find(abs(t_i - start_efit_time1) < tp/2);
    b = [];
    for i = 16:-1:1
        b = [b; sign_bb(16:-1:1, i, ind(1))];
    end
    b = double(b);
    b_data = [b_data, b];
end

% вычисляем верхнюю и нижнюю границе столбца свободных членов
N = length(b_data);
inf_b = zeros(N, 1); 
sup_b = zeros(N, 1); 
for i = 1:N
    inf_b(i) = min(b_data(i, :));
    sup_b(i)  = max(b_data(i, :));
end
%%
b = (sup_b + inf_b)/2;
A = hord_matrix;

%x = (A'A) * A' * b
disp(strcat("cond(A) = ", num2str(cond(A))));
disp(strcat("cond(A'A) = ", num2str(cond(A'*A))));
x2 = inv((A'*A)) *A' *b; 
lambda = eig(A'*A);
figure()
hist(lambda)  
title("собственные значения матрицы А'A")
ylabel("numder")
xlabel("lambda")

%tolsolvty
b_sup_tmp = sup_b;
b_inf_tmp = inf_b;
A_inf = A;
A_sup = A;
[tolmax, argmax, envs, ccode] =  tolsolvty(A_inf, A_sup, b_inf_tmp,  b_sup_tmp);
disp(strcat("tolmax = ", num2str(tolmax)));
iter = 0;
disp(strcat("iteration = ", num2str(iter)));

figure()
hold on;
grid on;    
plot_A = A_sup * argmax;
x_axis = 1:length(b_sup_tmp);
plot(x_axis, plot_A');
plot(x_axis, b_inf_tmp')
plot(x_axis, b_sup_tmp');
legend("Ax", "inf b", "sup b")
ylabel("value")
xlabel("i")
title("tolsolvty with K=2");

if(tolmax < 0)
    iter = iter + 1;
    disp(strcat("iteration = ", num2str(iter)));
    shift = -tolmax;
    disp(strcat("delta b = ", num2str(shift)));
    b_sup_tmp = b_sup_tmp + shift;
    b_inf_tmp = b_inf_tmp - shift;
    [tolmax, argmax, envs, ccode] =  tolsolvty(A_inf, A_sup, b_inf_tmp,  b_sup_tmp);
    disp(strcat("tolmax = ", num2str(tolmax)));

    figure()
    hold on;
    grid on;    
    % т.к. сейчас А - не интервальная матрица, то (A_inf == A_sup)
    % и для вычисления значения можем взять любую
    plot_A = A_sup * argmax;
    x_axis = 1:length(b_sup_tmp);
    plot(x_axis, plot_A');
    plot(x_axis, b_inf_tmp')
    plot(x_axis, b_sup_tmp');
    legend("Ax", "inf b", "sup b")
    ylabel("value")
    xlabel("i")
    title("tolsolvty with extended interval for b");
end


figure()
hold on;
grid on;
plot(argmax, "*");
ylabel("x_i")
xlabel("i")
title("tolsolvty decision")
figure()
title("tolsolvty histogram")
hold on;
grid on;
hist(argmax);
ylabel("value")
xlabel("x_i")

%%
%оценка числа обусловленности матрицы А
format short
matrix_radius = 0.1;
A_inf_1 = A * (1 -  matrix_radius);
A_sup_1 = A * (1 +  matrix_radius);
b_sup_tmp_1 = sup_b;
b_inf_tmp_1 = inf_b;

for i = 10:10:100
    cond_A = my_HeurMinCond(A_inf_1, A_sup_1, i);  
    disp(strcat("rad = ", num2str(matrix_radius)," : HeurMinCond(A, ", num2str(i), ") = ", num2str(cond_A)));
end
i = 1000;
cond_A = my_HeurMinCond(A_inf_1, A_sup_1, i);  
    disp(strcat("rad = ", num2str(matrix_radius)," : HeurMinCond(A, ", num2str(i), ") = ", num2str(cond_A)));
iter = 100;
for rad = 0.1:0.05:0.5
    A_inf_1 = A * (1 -  matrix_radius);
    A_sup_1 = A * (1 +  matrix_radius);
    b_sup_tmp_1 = sup_b;
    b_inf_tmp_1 = inf_b;
    cond_A = my_HeurMinCond(A_inf_1, A_sup_1, iter);        
    disp(strcat("rad = ", num2str(rad), " : HeurMinCond(A, ", num2str(iter), ") = ", num2str(cond_A)));
end

%%
%вычисление IVE интервальной матрицы А

%радиус элементов матрицы А - 10% от их величины
A_inf_1 = A * 0.9;
A_sup_1 = A * 1.1;
b_sup_tmp_1 = sup_b;
b_inf_tmp_1 = inf_b;

cond_A = my_HeurMinCond(A_inf_1, A_sup_1);

[tolmax_1, argmax_1, envs_1, ccode_1] =  tolsolvty(A_inf_1, A_sup_1, b_inf_tmp_1,  b_sup_tmp_1);
if(tolmax_1 < 0)
    shift_1 = -tolmax_1;
    b_sup_tmp_1 = b_sup_tmp_1 + shift_1;
    b_inf_tmp_1 = b_inf_tmp_1 - shift_1;
    [tolmax_1, argmax_1, envs_1, ccode_1] =  tolsolvty(A_inf_1, A_sup_1, b_inf_tmp_1,  b_sup_tmp_1);
end

A_IVE_1 = my_IVE(A_inf_1, A_sup_1,  b_inf_tmp_1,  b_sup_tmp_1, tolmax_1, argmax_1, length(argmax_1));
