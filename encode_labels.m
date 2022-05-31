function labels = encode_labels(xy_values, modulation_order, modulation_symbols)
    [~, n] = size(xy_values);
    labels = zeros(modulation_order, n);
    for i = 1:n
        xy_value =  xy_values(:, i);
        d = vecnorm((xy_value' - modulation_symbols)');
        [~, I] = min(d);
        label = -1 * ones(1, modulation_order);
        label(I) = 1;
        labels(:, i) = label;
    end
end