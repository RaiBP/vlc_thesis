function xy_values = rgb_signal_to_xy_values(rgb_signal, rgb_values)
    signal_length = length(rgb_signal);
    xy_values = zeros(2, signal_length);
    for i=1:signal_length
        Pi = rgb_signal(1, i);
        Pj = rgb_signal(2, i);
        Pk = rgb_signal(3, i);

        if Pi < 0 || Pj < 0 || Pk < 0
            Pmin = abs(min([Pi, Pj, Pk]));
            Pi = Pi + Pmin;
            Pj = Pj + Pmin;
            Pk = Pk + Pmin;
        end

        Pt = Pi+Pj+Pk;
        u = [Pi/Pt; Pj/Pt; Pk/Pt];  
        v = rgb_values'*u;
        xy_values(:, i) = [v(1); v(2)];
    end
end