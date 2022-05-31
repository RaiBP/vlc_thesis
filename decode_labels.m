function xy_values = decode_labels(labels, modulation_symbols, use_onehot_encoding)
    if use_onehot_encoding
        [~, numeric_labels] = max(labels);
    else
        numeric_labels = labels;
    end
    [~, n] = size(numeric_labels);
    xy_values = zeros(2, n);
    for i = 1:n
        xy_value = modulation_symbols(numeric_labels(i), :);
        xy_values(:, i) = xy_value;
    end
end