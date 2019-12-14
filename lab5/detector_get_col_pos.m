function [col_pos] = detector_get_col_pos(det, col)
POINT = point_get_interface();
p = POINT.create(0, 0);
little_step_num = floor(col/2);
big_step_num = floor((col - 1)/2);

%при равномерном распределении столбцов
%x_step = [det.horizontal_step.x, det.horizontal_step.x]
%y_step = [det.horizontal_step.y, det.horizontal_step.y]
%при неравномерном распределении столбцов
alpha = atan(det.direction.y / det.direction.x);
x_step = [det.step(1) * cos(alpha), det.step(2) * cos(alpha)];
y_step = [det.step(1) * sin(alpha), det.step(2) * sin(alpha)] ;

x_shift = little_step_num * x_step(1) + big_step_num * x_step(2);
y_shift = little_step_num * y_step(1) + big_step_num * y_step(2);
x = det.start.x + x_shift;
y = det.start.y + y_shift;
col_pos = POINT.create(x, y);
end