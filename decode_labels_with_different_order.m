function xy_values = decode_labels_with_different_order(labels, modulation_order, modulation_symbols, use_onehotencoding, label_order)
    if use_onehotencoding
        [~, numeric_label] = max(labels);
    else
        numeric_label = labels;
    end
    
    for i=1:modulation_order
        numeric_label(numeric_label == i) = -1*label_order(i); % we replace the original order by the order given in 'label_order'
    end
    numeric_label = abs(numeric_label);
    
    xy_values = decode_labels(numeric_label', modulation_symbols, 0);
    end