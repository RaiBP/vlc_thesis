function frames = form_frames(data, data_bits_per_frame, channel_estimation_sequence)
    frame_number = length(data) / data_bits_per_frame;
    frame_length = data_bits_per_frame + length(channel_estimation_sequence);
    frames = zeros(frame_number, frame_length);

    for i=1:frame_number
        data_in_frame =  data(1+(i-1)*data_bits_per_frame:i*data_bits_per_frame);
        frames(i, :) = [channel_estimation_sequence data_in_frame];
    end
end