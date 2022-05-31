function centroids = calculate_centroids(data, labels, number_of_clusters)
    [dim, ~] = size(data);
    centroids = zeros(dim, number_of_clusters);
    for i=1:number_of_clusters
        for l=1:dim
            centroids(l, i) = mean(data(l, labels==i));
        end
    end
end