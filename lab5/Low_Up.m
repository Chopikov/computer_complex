function [low,up, n_res] = Low_Up(start_efit_time,end_efit_time, t_i, tp, t_start, sign_bb)


efit_time_i = t_i(start_efit_time <= t_i & t_i <= end_efit_time);
ind_i = int32((efit_time_i - t_start) / tp + 1);
B_i = double(sign_bb(:, :, ind_i));
B_i = rot90(B_i, 2);
[cnt_rows, cnt_cols, n] = size(B_i);
B_low_i = B_i(:, 1:cnt_cols / 2, :);
B_up_i = B_i(:, cnt_cols / 2 + 1:cnt_cols, :);
sum_B_low = sum(sum(B_low_i));
sum_B_low = sum_B_low(:)';
sum_B_up = sum(sum(B_up_i));
sum_B_up = sum_B_up(:)';


low = sum_B_low;
up = sum_B_up;
n_res = n;
end

