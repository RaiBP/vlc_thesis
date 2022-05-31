function interframe_period_symbols = calculate_interframe_period_length(data_rate, bits_per_symbol, optical_clock_rate, interframe_period, padding_length)
    interframe_period_secs = interframe_period / optical_clock_rate;
    interframe_period_bits = data_rate * interframe_period_secs - padding_length;
    interframe_period_symbols = ceil(interframe_period_bits / bits_per_symbol);
end