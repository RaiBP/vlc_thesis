function xy_value = bits_to_xy_space(bits, modulation_order, bits_per_symbol, modulation_symbols)
    bit_sequences = int2bit(0:modulation_order-1, bits_per_symbol)';
    row_index = find(ismember(bit_sequences, bits', 'rows'));
    xy_value = modulation_symbols(row_index, :);
end