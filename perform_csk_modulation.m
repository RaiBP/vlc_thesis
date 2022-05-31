function modulated_frame = perform_csk_modulation(frame, modulation_order, bits_per_symbol, rgb_values, modulation_symbols)
    y = reshape(frame, bits_per_symbol, []);
    [~, L] = size(y);
    Pi_values = zeros(1, L);
    Pj_values = zeros(1, L);
    Pk_values = zeros(1, L);
    for i=1:L
        xy_value = bits_to_xy_space(y(:, i), modulation_order, bits_per_symbol, modulation_symbols);
        [Pi, Pj, Pk] = xy_values_to_rgb_power(xy_value, rgb_values);
        Pi_values(i) = Pi;
        Pj_values(i) = Pj;
        Pk_values(i) = Pk;
    end
    modulated_frame = [Pi_values; Pj_values; Pk_values];
end