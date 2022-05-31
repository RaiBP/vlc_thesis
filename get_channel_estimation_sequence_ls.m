function channel_estimation_sequence = get_channel_estimation_sequence_ls(modulation_order, modulation_symbols, rgb_values)
    % per IEEE standard
    walsh_code_i = [1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    walsh_code_j = [0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
    walsh_code_k = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1];
    rgb_signal = [walsh_code_i; walsh_code_j; walsh_code_k];
    xy_values = rgb_signal_to_xy_values(rgb_signal, rgb_values);
    channel_estimation_sequence = xy_values_to_bits(xy_values, modulation_order, log2(modulation_order), modulation_symbols);
end  