function DC_component = calculate_DC_component(channel, fs)
    spectrum = abs(fftshift(fft(channel)));
    precision = fs/length(channel);
    f = linspace(-fs/2+precision/2, fs/2-precision/2, length(channel));
    DC_component = spectrum(f==0);
end