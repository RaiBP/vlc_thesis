clear all;
clc;
close all;

number_of_samples = 1; % number of samples per symbol
data_rate = [24, 72, 96] * 10^6; % data rate of system for 4, 8 and 16-CSK respectively
Ptx = 1; % transmission power, in W

iterations = 10; % number of Monte Carlo iterations

EbN0=0:2:46; % SNR per bit vector
SNR = 10.^(EbN0./10); % SNR vector

rx_position_list = [3, 1, 1.8;
                    3, 1.5, 1.8;
                    4, 2, 1.8;
                    2.2, 2.5, 1.8;
                    1, 2.5, 1.6];

rx_position_index = 1;

tx_position = [3, 0.5, 4.5]; % TX position in (x, y, z) coordinates
rx_position = rx_position_list(rx_position_index, :); % RX in (x, y, z) coordinates

sig_length = 1e6;
max_PDSU_length = 65535;
data_bits_per_frame = (max_PDSU_length - 6) * 8; % MHR = 5 bytes, payload header = 1 byte, in total we remove 6 bytes
number_of_frames = ceil(sig_length / data_bits_per_frame);
sig_length = number_of_frames * data_bits_per_frame;

modulation_type = {'4csk', '8csk', '16csk'};
number_of_modulation_types = length(modulation_type);

equalization_type = {'ls', 'selm', 'sselm', 'uselm', 'ls+selm', 'ls+sselm', 'ls+uselm'};
number_of_equalization_types = length(equalization_type);

% [110-010-001] xy coordinates values as given in
% Yokoi, A., Son, l., & Bae, T. (2010). More description about CSK constellation. IEEE, 802(7), 3-46.
symbols_4csk = [0.19, 0.78;
                0.337, 0.393;
                0.09, 0.13;
                0.73, 0.27];

symbols_8csk = [0.19, 0.78;
                0.157, 0.563;
                0.37, 0.61;
                0.509, 0.396;
                0.09, 0.13;
                0.189, 0.326
                0.41, 0.2;
                0.73, 0.27];

symbols_16csk = [0.19, 0.78;
                0.239, 0.651;
                0.157, 0.563;
                0.37, 0.61;
                0.206, 0.434;
                0.337, 0.393;
                0.419, 0.481;
                0.386, 0.264;
                0.123, 0.347;
                0.172, 0.218;
                0.09, 0.13;
                0.303, 0.177;
                0.55, 0.44;
                0.599, 0.311;
                0.517, 0.223;
                0.73, 0.27];

% 
rgb_values = [0.73, 0.27; 0.19, 0.78; 0.09, 0.13];

% from Chris R. Schoenenman (1991), Converting Frequency to RGB, comp.graphics.algorithms. http://steve.hollasch.net/cgindex/color/freq-rgb.html
xy_spectrum = [0.1741, 0.0050; 0.1740, 0.0050; 0.1738, 0.0049; 0.1736, 0.0049;
            0.1733, 0.0048; 0.1730, 0.0048; 0.1726, 0.0048; 0.1721, 0.0048;
            0.1714, 0.0051; 0.1703, 0.0058; 0.1689, 0.0069; 0.1669, 0.0086;
            0.1644, 0.0109; 0.1611, 0.0138; 0.1566, 0.0177; 0.1510, 0.0227;
            0.1440, 0.0297; 0.1355, 0.0399; 0.1241, 0.0578; 0.1096, 0.0868;
            0.0913, 0.1327; 0.0687, 0.2007; 0.0454, 0.2950; 0.0235, 0.4127;
            0.0082, 0.5384; 0.0039, 0.6548; 0.0139, 0.7502; 0.0389, 0.8120;
            0.0743, 0.8338; 0.1142, 0.8262; 0.1547, 0.8059; 0.1929, 0.7816;
            0.2296, 0.7543; 0.2658, 0.7243; 0.3016, 0.6923; 0.3373, 0.6589;
            0.3731, 0.6245; 0.4087, 0.5896; 0.4441, 0.5547; 0.4788, 0.5202;
            0.5125, 0.4866; 0.5448, 0.4544; 0.5752, 0.4242; 0.6029, 0.3965;
            0.6270, 0.3725; 0.6482, 0.3514; 0.6658, 0.3340; 0.6801, 0.3197;
            0.6915, 0.3083; 0.7006, 0.2993; 0.7079, 0.2920; 0.7140, 0.2859;
            0.7190, 0.2809; 0.7230, 0.2770; 0.7260, 0.2740; 0.7283, 0.2717;
            0.7300, 0.2700; 0.7311, 0.2689; 0.7320, 0.2680; 0.7327, 0.2673;
            0.7334, 0.2666; 0.7340, 0.2660; 0.7344, 0.2656; 0.7346, 0.2654;
            0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653;
            0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653;
            0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653;
            0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653; 0.7347, 0.2653;
            0.7347, 0.2653];

wavelengths = 380:5:780;

closest_wavelengths = zeros(1, 3);
for i=1:3
    closest_wavelengths(i) = wavelengths(find(min(vecnorm((xy_spectrum-rgb_values(i, :))')) == vecnorm((xy_spectrum-rgb_values(i, :))')));
end

pd_responsivities = get_photodiode_responsivities(closest_wavelengths);

% LEEFilters: https://leefilters.com/lighting/colour-effect-lighting-filters/
% Used in E. Monteiro & S. Hranilovic (2014). Design and implementation of color-shift keying for visible light communications. Journal of Lightwave Technology, 32(10), pp. 2053-2060.
transmission_spectra = [80.351, 1.428, 1.908;   % R: 24 Scarlet
                        11.769, 75.148, 13.076; % G: 738 JAS Green
                        0.279, 27.244, 70.684] / 100; % B: 141 Bright Blue
                    %   R       G       B

net_responsivities = transmission_spectra .* pd_responsivities';

modulation_symbols = {symbols_4csk, symbols_8csk, symbols_16csk};

for i=1:length(modulation_symbols)
    [modulation_order, ~] = size(modulation_symbols{i});
    channel_estimation_sequence_elm{i} = get_channel_estimation_sequence_elm(modulation_order, modulation_symbols{i});
    channel_estimation_sequence_uselm{i} = xy_values_to_bits(modulation_symbols{i}', modulation_order, log2(modulation_order), modulation_symbols{i});
    channel_estimation_sequence_ls{i} = get_channel_estimation_sequence_ls(modulation_order, modulation_symbols{i}, rgb_values);

end

optical_clock = 24e6;
interframe_period = 400; % LIFS = 400, SIFS = 120, RIFS = 40

equalization_parameters = get_equalization_parameters(modulation_type, equalization_type);
SNR_length = length(SNR);
for iter=1:iterations
    fprintf("Iteration %d\n", iter);
    [h, fs_CIR] = get_channel_model(tx_position, rx_position);
    H0 = calculate_DC_component(h, fs_CIR);
    delay = round(mean(grpdelay(h)));

    ber_list = zeros(number_of_equalization_types, number_of_modulation_types, SNR_length);
    

    for i=1:SNR_length
        tic
        fprintf("SNR: %d dB\n", 10*log10(SNR(i)));

        N0xB = (net_responsivities*H0*Ptx)^2 / SNR(i); % total noise power

        sgma = sqrt(N0xB); % standard deviation of noise after matched filter

        data = round(rand(1, sig_length));
        
        for k=1:number_of_equalization_types
          
            switch equalization_type{k}
                case 'ls'
                    channel_estimation_sequences = channel_estimation_sequence_ls; % LS channel estimation sequence
                case 'uselm'
                    channel_estimation_sequences = channel_estimation_sequence_uselm; % USELM channel estimation sequence
                otherwise
                    channel_estimation_sequences = channel_estimation_sequence_elm; % ELM channel estimation sequence
            end
            for l=1:number_of_modulation_types
                modulation_order = length(modulation_symbols{l});
                bits_per_symbol = log2(modulation_order);
                channel_estimation_sequence_modulated = perform_csk_modulation(channel_estimation_sequences{l}, modulation_order, bits_per_symbol, rgb_values, modulation_symbols{l});
                frames = form_frames(data, data_bits_per_frame, channel_estimation_sequences{l});
                [modulated_signal, indices] = modulate_frames(frames, length(channel_estimation_sequences{l}), modulation_symbols{l}, data_rate(l), number_of_samples, optical_clock, rgb_values, interframe_period);
                received_signal = apply_channel(modulated_signal, h, net_responsivities, sgma, data_rate(l), fs_CIR, number_of_samples, delay, modulation_order);
                [frame_matrix_rx, data_matrix_rx, estimation_sequence_matrix_rx] = extract_frames(received_signal, indices, number_of_frames);
                equalized_data = equalization(data_matrix_rx, estimation_sequence_matrix_rx, channel_estimation_sequence_modulated, number_of_frames, number_of_samples, modulation_order, equalization_type{k}, bits_per_symbol, modulation_symbols{l}, equalization_parameters{k, l}, data_bits_per_frame, rgb_values);
                [~, ber] = biterr(equalized_data, data);
                ber_list(k, l, i) = ber;
            end
          
        end
        toc
    end
    filename = sprintf('ber_iter%d_position%d.mat', iter, rx_position_index);
    save(filename, 'ber_list');
end

