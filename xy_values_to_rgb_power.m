function [Pi, Pj, Pk] = xy_values_to_rgb_power(xy_value, rgb_values)
    A = [rgb_values'; 1, 1, 1];
    u = [xy_value'; 1];
    v = pinv(A)*u;
    Pi = v(1);
    Pj = v(2);
    Pk = v(3);
end