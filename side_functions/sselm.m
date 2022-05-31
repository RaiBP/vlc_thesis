function elmModel=sselm(Xl,Yl,Xu,L,paras)

[l,elmModel.InputDim]=size(Xl);
u=size(Xu,1);
N=l+u;


% Decide whether it is a binary or multi-class problem
elmModel.labs=unique(Yl);
if length(elmModel.labs)==2     % Binary classification
    elmModel.MultiOutput=0;
    elmModel.OutputDim=1;
else                      % Multi-class classification
    elmModel.MultiOutput=1;
    elmModel.OutputDim=length(elmModel.labs);
end

% Random generate input weights
elmModel.InputWeight=rand(elmModel.InputDim,paras.NumHiddenNeuron)*2-1;

% Calculate hidden neuron output matrix
elmModel.Kernel=paras.Kernel;
switch paras.Kernel
    case 'sigmoid'
        H=1 ./ (1 + exp(-[Xl;Xu]*elmModel.InputWeight));
    case 'rbf'
        % to be added
end
Hl=H(1:l,:);
clear Xl Xu


% Calculate output weights
if elmModel.MultiOutput==0  % Binary classification
    Y=[Yl;zeros(u,1)];
    Cl_diag=l/2*((Yl==max(Yl))/sum(Yl==max(Yl))+(Yl==min(Yl))/sum(Yl==min(Yl)));
else                                % Multi-class classification
    Y=zeros(N,elmModel.OutputDim);
    Y(1:l,:)=-1;    
    Cl_diag=zeros(l,1);
    for i=1:elmModel.OutputDim
        Y(Yl==elmModel.labs(i),i)=1;
        Cl_diag=Cl_diag+(Yl==elmModel.labs(i))/(sum(Yl==elmModel.labs(i))*elmModel.OutputDim/l);
    end
end

Cl=diag(paras.C*Cl_diag);
if  (paras.NumHiddenNeuron>N) % Comput C only the training instances is less than hidden nodes
    C=diag([paras.C*Cl_diag;zeros(u,1)]); 
end


t_elm_start=tic;
if  (paras.NumHiddenNeuron<N)
    elmModel.OutputWeight=(eye(paras.NumHiddenNeuron)+Hl'*Cl*Hl+paras.lambda*H'*L*H)\ (Hl'*Cl* Y(1:l,:));
else
    A=eye(N)+(C+paras.lambda*L)*(H*H');
    B=C*Y;
    D=A\B;
    elmModel.OutputWeight=H'*D;
    elmModel.OutputWeight0=H'*((eye(N)+(C+paras.lambda*L)*(H*H'))\(C*Y));
    norm(elmModel.OutputWeight-elmModel.OutputWeight0)
end
elmModel.TrainTime=toc(t_elm_start);

% Calculate training accuracy
if elmModel.MultiOutput==0  % Binary classification
    TrainAccuracy=100*mean(sign(Hl*elmModel.OutputWeight)==Yl);
else                                % Multi-class classification
    [~,idx]=max((Hl*elmModel.OutputWeight)');
    TrainAccuracy=100*mean(elmModel.labs(idx)==Yl);
end

if ~paras.NoDisplay
    disp(['Traning time is ',num2str(elmModel.TrainTime)])
    disp(['Traning accracy is ',num2str(TrainAccuracy),'%'])
end


