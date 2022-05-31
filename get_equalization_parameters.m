function equalization_parameters = get_equalization_parameters(modulation_type, equalization_type)
    for i=1:length(equalization_type)
        eq_type = equalization_type{i};
        if strcmp(eq_type, 'ls') == 0
            if contains(eq_type, 'ls+')
                eq_type = extractAfter(eq_type, 'ls+');
            end
            for l=1:length(modulation_type)
                filename = sprintf("gridsearch_%s_25dB_%s.csv", eq_type, modulation_type{l});
                filepath = fullfile('.','equalization_params', filename);
                params_matrix = readmatrix(filepath);
                [~, cols_number] = size(params_matrix);
                params_matrix_sorted = sortrows(params_matrix, cols_number); % sort by BER (last column)
                for k=1:cols_number-1
                    params{k} = params_matrix_sorted(1, k);
                end
                equalization_parameters{i, l} = params;
            end
        end
    end
end