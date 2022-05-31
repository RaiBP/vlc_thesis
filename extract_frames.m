function [frame_matrix, data_matrix, estimation_sequence_matrix] = extract_frames(signal, indices, number_of_frames)
    for i=1:number_of_frames
        start_index = indices(1, i);
        end_index = indices(2, i);
        channel_index = indices(3, i);

        frame_matrix{i} = signal(:, start_index:end_index);
        data_matrix{i} = signal(:, channel_index+1:end_index);
        estimation_sequence_matrix{i} = signal(:, start_index:channel_index);
    end
end