function x_output = selm_equalization(y_train, x_train, x_test, NN)
NumberofHiddenNeurons=NN;

size_train = size(x_train);
size_test = size(x_test);

NumberofTrainingData=size_train(2);
NumberofTestingData=size_test(2);
NumberofInputNeurons=size_train(1);
InputWeight=rand(NumberofHiddenNeurons,NumberofInputNeurons)*2-1;

BiasofHiddenNeurons=rand(NumberofHiddenNeurons,1);
tempH=InputWeight*x_train+BiasofHiddenNeurons(:,ones(1,NumberofTrainingData));

H = 1 ./ (1 + exp(-tempH));

OW=pinv(H') * y_train';

tempH_test=InputWeight*x_test + BiasofHiddenNeurons(:,ones(1,NumberofTestingData));
H_test = 1 ./ (1 + exp(-tempH_test));
x_output=(H_test' * OW)';
end