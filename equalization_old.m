function equilized_data = equalization_old(data_matrix, estimation_sequences_rx, estimation_sequences_tx, modulation_order, equalization_type, bits_per_symbol, modulation_symbols, parameters, rgb_values)
    equilized_data = [];
    estimation_sequence_tx_xy = 0;
    if strcmp(equalization_type, "ls") == 0
        estimation_sequence_tx_xy = rgb_signal_to_xy_values(real_estimation_sequence, rgb_values);
    end

    for i=1:number_of_frames
        received_estimation_sequence = estimation_sequences_rx(i, :);
        received_data = data_matrix(i, :);

        received_estimation_sequence_downsamp = mean(reshape(received_estimation_sequence, number_of_samples, []));
        received_data_downsamp = mean(reshape(received_data, number_of_samples, []));

        if strcmp(equalization_type, "ls") > 0
            equilized_data_bits = ls_equalization(estimation_sequences_tx, received_estimation_sequence_downsamp, received_data_downsamp, parameters);
            
        else
            received_data_xy = rgb_signal_to_xy_values(received_data_downsamp, rgb_values);
            estimation_sequence_rx_xy = rgb_signal_to_xy_values(received_estimation_sequence_downsamp, rgb_values);
            
            switch equalization_type
                case "relm"
                    equalized_data_xy = relm_equalization(estimation_sequence_tx_xy, estimation_sequence_rx_xy, received_data_xy, parameters);
                case "sselm"
                    equalized_data_xy = sselm_equalization(estimation_sequence_tx_xy, estimation_sequence_rx_xy, received_data_xy, parameters);
                case "uselm"
                    equalized_data_xy = uselm_equalization(estimation_sequence_tx_xy, estimation_sequence_rx_xy, received_data_xy, parameters);
                otherwise
                    equalized_data_xy = elm_equalization(estimation_sequence_tx_xy, estimation_sequence_rx_xy, received_data_xy, parameters);
            end
            equilized_data_bits = xy_values_to_bits(equalized_data_xy, modulation_order, bits_per_symbol, modulation_symbols);
        end    
        equilized_data = [equilized_data equilized_data_bits];
    end
    
end