function [acc,MSE,auc,predict]=sselm_predict(X,Y,elmModel)


% Calculate hidden neuron output matrix
switch elmModel.Kernel
    case 'sigmoid'
        H=1 ./ (1 + exp(-X*elmModel.InputWeight));
    case 'rbf'
        H=calckernel(elmModel,elmModel.InputWeight,X);
end


% Calculate training accuracy
if elmModel.MultiOutput==0  % Binary classification
    out=H*elmModel.OutputWeight;
    acc=100*mean(sign(H*elmModel.OutputWeight)==Y);
    [~,~,~,auc] = perfcurve(Y,out,1);
    predict=sign(out);
    MSE=mse(Y-out);
else                                % Multi-class classification
    y=-ones(length(Y),elmModel.OutputDim);
    for i=1:elmModel.OutputDim
        y(Y==elmModel.labs(i),i)=1;
    end
    out=H*elmModel.OutputWeight;
    [~,idx]=max(out');
    predict=elmModel.labs(idx);
    acc=100*mean(predict==Y);
    auc=NaN;
    MSE=mse(y-out);
end




