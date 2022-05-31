function x_output_labels = sselm_equalization(pilot_tx, pilot_rx, frame_rx, NumberofHiddenNeurons, C, lambda, lu, pilot_length)
    addpath(genpath('side_functions'))
    options.NN=7;%17qpsk,7qam
    options.GraphWeights='binary';
    options.GraphDistanceFunction='euclidean';
    options.LaplacianNormalize=1;
    options.LaplacianDegree=1;
    paras.NumHiddenNeuron=NumberofHiddenNeurons;
    paras.C=C;
    paras.lambda=lambda;
    paras.Kernel='sigmoid';
    paras.NoDisplay=0;
    T=pilot_tx(:,1:pilot_length)';
    P=pilot_rx(:,1:pilot_length)';
    PP=frame_rx';
    Xu = PP(1:lu, :);
    L=laplacian(options,[P; Xu]);
    elmModel=sselm(P,T,Xu,L,paras);
    %Htra=1 ./ (1 + exp(-P*elmModel.InputWeight));
    %x_output_wp=Htra*elmModel.OutputWeight;
    H=1 ./ (1 + exp(-frame_rx'*elmModel.InputWeight));
    x_output=H*elmModel.OutputWeight;
    [~,x_output_labels]=max(x_output');
end