function received_signal = apply_channel(signal, channel, responsitivy_matrix, noise, data_rate, CIR_sampling_frequency, number_of_samples, delay, modulation_order)
    bits_per_symbol = log2(modulation_order);
    signal_sampling_frequency = number_of_samples * data_rate / bits_per_symbol;
    upsampling_factor = floor(CIR_sampling_frequency / signal_sampling_frequency);
    upsampled_number_of_samples = number_of_samples * upsampling_factor;

    for i=1:3
        signal_single_color = signal(i, :);
        signal_length = length(signal_single_color);
        indeces = 1:1/upsampled_number_of_samples:signal_length;
        signal_interp = interp1(1:signal_length, signal_single_color, indeces); % linear interpolation
        rx_signal_interp = conv(channel, signal_interp);
        rx_signal = rx_signal_interp(1+delay:upsampled_number_of_samples:end);
        rx_signal_cell{i} = max(rx_signal, 0);
    end

    for i=1:3
        rx_signal_multichannel{i} = 0;
        for l=1:3
            rx_signal_multichannel{i} = rx_signal_multichannel{i} + rx_signal_cell{l}*responsitivy_matrix(i, l)*2 + noise(i, l)*randn(1, length(rx_signal_cell{l}));
        end
        rx_signal_multichannel{i} = max(rx_signal_multichannel{i}, 0);
    end

    received_signal = [rx_signal_multichannel{1}; rx_signal_multichannel{2}; rx_signal_multichannel{3}];
end