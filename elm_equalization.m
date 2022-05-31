function equilized_data = elm_equalization(data_matrix, received_estimation_sequences, real_estimation_sequence, modulation_order, number_of_samples, parameters)
    number_of_neurons = parameters(1);
    estimation_sequence_tx_xy = signal_to_xy_values(real_estimation_sequence);
    
    equilized_data = [];
    for i=1:number_of_frames
        received_estimation_sequence = received_estimation_sequences(i, :);
        received_data = data_matrix(i, :);

        received_estimation_sequence_downsamp = mean(reshape(received_estimation_sequence, number_of_samples, []));
        received_data_downsamp = mean(reshape(received_data, number_of_samples, []));
        
        received_data_xy = signal_to_xy_values(received_data_downsamp);
        estimation_sequence_rx_xy = signal_to_xy_values(received_estimation_sequence_downsamp);
        
        received_data_elm_xy = elm_equalization(estimation_sequence_tx_xy, estimation_sequence_rx_xy, received_data_xy, number_of_neurons);

        received_data_elm_bits = xy_values_to_bits(received_data_elm_xy);
        equilized_data = [equilized_data received_data_elm_bits]
    end
end