function responsivity = get_photodiode_responsivities(wavelengths)
    % Data from Thorlabs PDA36A Si Switchable Gain Detector
    % Found in https://www.thorlabs.com/drawings/a24ce0c193c76ec0-187393EC-DE5C-068A-EB49EEF1A1159D8E/PDA36A-Manual.pdf
    % Used in E. Monteiro & S. Hranilovic (2014). Design and implementation of color-shift keying for visible light communications. Journal of Lightwave Technology, 32(10), pp. 2053-2060.
    responsivity_data = readmatrix('photodiode_responsivity.csv', 'Delimiter', ';', 'DecimalSeparator', ',');
    wavelength_axis = responsivity_data(:, 1);
    responsivity_axis = responsivity_data(:, 2);
    responsivity = zeros(1, 3);
    for i=1:3
        responsivity(i) = responsivity_axis(find(min(abs(wavelength_axis-wavelengths(i))) == abs(wavelength_axis-wavelengths(i))));
    end
end