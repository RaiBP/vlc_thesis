function equalized_data = ls_equalization(estimation_sequence_tx, estimation_sequence_rx, data_rx)
    estimation_sequence_rx_averaged = squeeze(mean(reshape(estimation_sequence_rx, 3, 2, []), 2));
    estimation_sequence_tx_averaged = squeeze(mean(reshape(estimation_sequence_tx, 3, 2, []), 2));
    channel_estimation = estimation_sequence_rx_averaged * pinv(estimation_sequence_tx_averaged);
    equalized_data = inv(channel_estimation) * data_rx;
end