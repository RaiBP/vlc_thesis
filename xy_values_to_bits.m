function bits = xy_values_to_bits(xy_values, modulation_order, bits_per_symbol, symbols)
    bits = [];
    for i=1:length(xy_values)
        d = vecnorm((xy_values(:, i)' - symbols)');
        [~, I] = min(d);
        bit_sequences = int2bit(0:modulation_order-1, bits_per_symbol)';
        bits = [bits bit_sequences(I, :)];
    end
end