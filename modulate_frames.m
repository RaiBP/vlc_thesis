function [upsampled_signal, indices] = modulate_frames(frames, estimation_sequence_length, modulation_symbols, data_rate, number_of_samples, optical_clock_rate, rgb_values, interframe_period)
    M = length(modulation_symbols);
    m = log2(M);
    [frame_number, frame_length] = size(frames);
    padding_length = m - mod(frame_length, m);
    padding = zeros(1, padding_length);

    frame_start_indices = zeros(1, frame_number);
    frame_end_indices = zeros(1, frame_number);
    estimation_sequence_indices = zeros(1, frame_number);

    frame_start_index = 1;

    modulated_signal = [];
    for i=1:frame_number

        frame = [frames(i, :) padding];
        frame_end_index = frame_start_index + ceil((length(frame))/m);
        estimation_sequence_index = frame_start_index + ceil(estimation_sequence_length/m);

        frame_start_indices(i) = frame_start_index;
        frame_end_indices(i) = frame_end_index;
        estimation_sequence_indices(i) = estimation_sequence_index;

        modulated_frame = perform_csk_modulation(frame, M, m, rgb_values, modulation_symbols);
        interframe_length_symbols = calculate_interframe_period_length(data_rate, m, optical_clock_rate, interframe_period, padding_length);
        interframe_period_modulated = zeros(3, interframe_length_symbols);
        modulated_signal = [modulated_signal modulated_frame interframe_period_modulated];
        
        frame_start_index = frame_end_index + interframe_length_symbols;
    end

    upsampled_signal = zeros(3, number_of_samples*length(modulated_signal));
    for l=1:3
        upsampled_signal(l, :) = rectpulse(modulated_signal(l, :), number_of_samples);
    end

    indices = [number_of_samples*(frame_start_indices-1)+1; number_of_samples*(frame_end_indices-1); number_of_samples*(estimation_sequence_indices-1)];
end