function [x_lim, y_lim]=get_plot_lim(x, y)

min_x_sep = min(x);
max_x_sep = max(x);
dx = max_x_sep - min_x_sep;
min_y_sep = min(y);
max_y_sep = max(y);
dy = max_y_sep - min_y_sep;
if(dx < dy)
    half = (max_y_sep - min_y_sep) / (2 * 0.9); %чтобы график не прилипал к краю оставим зазор в 20%
else
    half = (max_x_sep - min_x_sep) / (2 * 0.9); %чтобы график не прилипал к краю оставим зазор в 20%
end
    mid_x = (min_x_sep + max_x_sep) / 2;
    x_lim_left = mid_x - half;
    x_lim_right = mid_x + half;
    mid_y = (min_y_sep + max_y_sep) / 2;
    y_lim_left = mid_y - half;
    y_lim_right = mid_y + half;
    
    x_lim = [x_lim_left, x_lim_right];
    y_lim = [y_lim_left, y_lim_right];
    
end