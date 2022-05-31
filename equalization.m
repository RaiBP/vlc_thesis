function equilized_data = equalization(data_matrix, estimation_sequences_rx, estimation_sequences_tx, number_of_frames, number_of_samples, modulation_order, equalization_type, bits_per_symbol, modulation_symbols, parameters, data_length, rgb_values)
    equilized_data = [];
    pre_equalization = 0;

    for i=1:number_of_frames
        received_estimation_sequence = estimation_sequences_rx{i};
        received_data = data_matrix{i};

        received_estimation_sequence_downsamp = squeeze(mean(reshape(received_estimation_sequence, 3, number_of_samples, []), 2));
        received_data_downsamp = squeeze(mean(reshape(received_data, 3, number_of_samples, []), 2));

        if strcmp(equalization_type, "ls") > 0
            equilized_data_rgb = ls_equalization(estimation_sequences_tx, received_estimation_sequence_downsamp, received_data_downsamp);
            equalized_data_xy = rgb_signal_to_xy_values(equilized_data_rgb, rgb_values);
         
        else
            if contains(equalization_type, "ls+") > 0
                equalization_type = extractAfter(equalization_type, "ls+");
                [received_data_downsamp, received_estimation_sequence_downsamp] = ls_preequalization(estimation_sequences_tx, received_estimation_sequence_downsamp, received_data_downsamp);
                pre_equalization = 1;
            end

            labels_tx = encode_labels(rgb_signal_to_xy_values(estimation_sequences_tx, rgb_values), modulation_order, modulation_symbols);
            use_onehot_encoding = 0;
            switch equalization_type
                case "relm"
                    use_onehot_encoding = 1;
                    number_of_neurons = parameters{1};
                    gamma = parameters{2};
                    equalized_data_labels = relm_equalization(labels_tx, received_estimation_sequence_downsamp, received_data_downsamp, number_of_neurons, gamma);
                case "sselm"
                    number_of_neurons = parameters{1};
                    C = parameters{2};
                    lambda = parameters{3};
                    [~, numeric_labels_tx] = max(labels_tx);
                    equalized_data_labels = sselm_equalization(numeric_labels_tx, received_estimation_sequence_downsamp, ...
                                                                received_data_downsamp, number_of_neurons, C, lambda, ...
                                                                round(length(received_data_downsamp)/25), length(labels_tx));
                case "uselm"
                    number_of_neurons = parameters{1};
                    embedding_dimension = parameters{2};
                    lambda = parameters{3};
                    equalized_data_labels = uselm_equalization(received_data_downsamp, number_of_neurons, lambda, embedding_dimension, ...
                                                                round(length(received_data_downsamp)/25), modulation_order);

                otherwise
                    use_onehot_encoding = 1;
                    number_of_neurons = parameters{1};
                    equalized_data_labels = selm_equalization(labels_tx, received_estimation_sequence_downsamp, received_data_downsamp, number_of_neurons);     
            end
            
            if strcmp(equalization_type, "uselm") > 0
                centroids_rgb = calculate_centroids(received_data_downsamp, equalized_data_labels, modulation_order);
                centroids_xy = rgb_signal_to_xy_values(centroids_rgb, rgb_values);
                if pre_equalization
                    centroid_labels = encode_labels(centroids_xy, modulation_order, modulation_symbols);
                else
                    modulation_symbols_rx_xy = rgb_signal_to_xy_values(received_estimation_sequence_downsamp, rgb_values);
                    centroid_labels = encode_labels(centroids_xy, modulation_order, modulation_symbols_rx_xy');
                end
                [~, numeric_centroid_labels] = max(centroid_labels);
                equalized_data_xy = decode_labels_with_different_order(equalized_data_labels, modulation_order, modulation_symbols, use_onehot_encoding, numeric_centroid_labels);

            else
                equalized_data_xy = decode_labels(equalized_data_labels, modulation_symbols, use_onehot_encoding);
            end
        end
        equilized_data_bits = xy_values_to_bits(equalized_data_xy, modulation_order, bits_per_symbol, modulation_symbols);
        equilized_data = [equilized_data equilized_data_bits(1:data_length)];
    end
    
end