function [equalized_data, equalized_channel_estimation] = ls_preequalization(estimation_sequence_tx, estimation_sequence_rx, data_rx)
    channel_estimation = estimation_sequence_rx * pinv(estimation_sequence_tx);
    channel_estimation_inv = inv(channel_estimation);
    equalized_data = channel_estimation_inv * data_rx;
    equalized_channel_estimation = channel_estimation_inv * estimation_sequence_rx;
end