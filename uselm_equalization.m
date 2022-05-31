function label_uselm = uselm_equalization(data, NN, lambda, embedding_dimensions, training_data_length, modulation_order)
    addpath(genpath('side_functions'))
    
    X = data'; % input
    X_train = data(:, 1:training_data_length)';
    
    NC=modulation_order; % specify number of clusters
    
    % %%%%%%%%%%%%%%%%% Step 1: construct graph Laplacian %%%%%%%%%%%%%%%%%
    % hyper-parameter settings for graph
    options.GraphWeights='binary';
    options.GraphDistanceFunction='euclidean';
    options.LaplacianNormalize=0;
    options.LaplacianDegree=1;
    options.NN=5;
    
    L=laplacian(options, X_train);
    
    paras.NE=embedding_dimensions; % specify dimensions of embedding
    paras.NumHiddenNeuron=NN;
    paras.NormalizeInput=0;
    paras.NormalizeOutput=0;
    paras.Kernel='sigmoid';
    paras.lambda=lambda;
    elmModel = uselm(X, X_train, L, paras);
        [label_uselm, ~] = litekmeans(elmModel.Embed, NC, 'MaxIter', 200);
    end