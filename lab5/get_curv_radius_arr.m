function r = get_curv_radius_arr(RBDRY, ZBDRY, NBDRY)
    %вычисляем отдельно для первой и последней точки, учитывая замкнутость кривой
    a = [RBDRY(NBDRY),ZBDRY(NBDRY)];
    b = [RBDRY(1),ZBDRY(1)];
    c = [RBDRY(2),ZBDRY(2)];      
    r = get_curv_radius(a, b, c);
    for i=2:NBDRY - 1
        a = [RBDRY(i-1),ZBDRY(i-1)];
        b = [RBDRY(i),ZBDRY(i)];
        c = [RBDRY(i + 1),ZBDRY(i + 1)];
        tmp_r = get_curv_radius(a, b, c);
        r = [r, tmp_r];
    end

    a = [RBDRY(NBDRY-1),ZBDRY(NBDRY-1)];
    b = [RBDRY(NBDRY),ZBDRY(NBDRY)];
    c = [RBDRY(1),ZBDRY(1)];
    tmp_r = get_curv_radius(a, b, c);
    r = [r, tmp_r];
end