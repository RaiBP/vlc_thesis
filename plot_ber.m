
EbN0=0:2:46;
% signal-to-noise ratio in dB.
SNR = 10.^(EbN0./10);

n_iters = 4;
modulation_type = {'4csk', '8csk', '16csk'};
equalization_type = {'ls', 'selm', 'sselm', 'uselm', 'ls+selm', 'ls+sselm', 'ls+uselm'};

% eq, mod, snr
sums = num2cell(zeros(length(equalization_type), length(modulation_type)));
for iter=1:n_iters
    temp = load(sprintf("ber_iter%d_position1.mat", iter));
    ber_list = temp.ber_list;
    for i=1:length(equalization_type)
        for k=1:length(modulation_type)
            sums{i, k} = sums{i, k} + ber_list(i, k, :);
        end
    end
end

k = 3;
for i=1:length(equalization_type)
        semilogy(EbN0, squeeze(sums{i, k})./n_iters, '-*','markersize',2,'linewidth',1,'DisplayName',sprintf('%s, %s', equalization_type{i}, modulation_type{k}));
        hold on
end
semilogy(EbN0,qfunc(sqrt(10.^(EbN0/10))),'k-x','markersize',4,'linewidth',2,'DisplayName','AWGN channel');
grid on
 legend show
 xlim([0 46])
 %ylim([1e-2 1])
 title('16-CSK')
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');