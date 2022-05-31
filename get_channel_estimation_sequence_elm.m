function channel_estimation_sequence = get_channel_estimation_sequence_elm(modulation_order, modulation_symbols)
    sequence_length_in_symbols = 24; % per IEEE standard
    symbol_repetitions = floor(sequence_length_in_symbols / modulation_order);
    remainder = mod(sequence_length_in_symbols, modulation_order);
    channel_estimation_sequence_xy_values = repmat(modulation_symbols', [1, symbol_repetitions]);
    if remainder ~= 0
        channel_estimation_sequence_xy_values = [channel_estimation_sequence_xy_values modulation_symbols(1:remainder, :)'];
    end
    channel_estimation_sequence = xy_values_to_bits(channel_estimation_sequence_xy_values, modulation_order, log2(modulation_order), modulation_symbols);
end  