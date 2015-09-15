% Warped stretch group delay design by spectrotemporal distribution
% By Ata Mahjoubfar and Claire Chen at Jalali-Lab, UCLA
clear all
clc
lambda = 1.5e-6; % Carrier wavelength [m]
c = 299792458; % Speed of light [m/s]
f0 = c/lambda; % Carrier frequency [Hz]

Fs = 1e13; % Sampling rate [Hz]
dt = 1/Fs; % Time simulation resolution [s]
max_t = 0.333e-8; % Maximum simulation time [s]
t = -max_t : dt : max_t; % Time samples [s]
N = length(t); % Number of samples
df = 0.5/max_t; % Frequency simulation resolution [Hz]
f = -0.5/dt : df : 0.5/dt; % Frequency samples [Hz]
f_GD = (f(2:end)+f(1:end-1))/2; % Frequencies of GD samples
f_GDD = f(2:end-1); % Frequencies of GDD samples
added_chirp = exp(+1j*1e-13*pi*f.^2*(max_t)) .* exp(-1j*2e-26*pi*f.^3*(max_t)); % Input signal chirp
asymmetry = true; % Input spectrum asymmetry
switch 'set temporal shape'
    case 'set temporal shape'
        E_in_t = 1.5e3.*(exp(-((t-7e-12)/1e-12).^2)+0.7*exp(-((t+7e-12)/1e-12).^2)+10*exp(-(t/2e-13).^2)); % Complex envelope of input electric field in time [V/m]
        E_in_f = 1/Fs * added_chirp.*fftshift(fft(E_in_t)); % Complex envelope of input electric field in frequency [V/m/Hz]
        E_in_t = Fs * ifft(ifftshift(E_in_f));
    case 'set spectral shape'
        E_in_f = added_chirp.*((0.25.*exp(-(f/0.5e12).^2) + ((0.5 + double(~asymmetry).*0.5*cos(2*pi*(1/150e9)*f)).*exp(-((f-(1.5e12 - double(asymmetry).*0.6e12))/0.15e12).^2) + (0.5 + 0.5*cos(2*pi*(1/150e9)*f)).*exp(-((f+1.5e12)/0.15e12).^2))).*1e-8.*exp(1j*2*pi*f*(max_t)));
        E_in_t = Fs * ifft(ifftshift(E_in_f));
end
E_in_envelope_t = abs(E_in_t); % Envelope of input electric field in time [V/m]
E_in_envelope_f = 1/Fs * fftshift(fft(E_in_envelope_t)); % Envelope of input electric field in frequency [V/m/Hz]
E_in_unchirped_t = Fs * ifft(ifftshift(exp(1j*2*pi*f*max_t).*abs(E_in_f))); % Unchirped input electric field in time [V/m]

c1 = 2*pi*600*(1e-12).^2;
c2 = 2*pi*2*(1e-12).^4;
c3 = 1e12;
c4 = c1*c3;
c5 = 1.4e12;
c7 = 1.3*c5;
c8 = -c1*c7^2*pi^2/8;
filter_phases_wihtout_chirp = {@(f) c1/2.*f.^2 ...
    @(f) c4.*(f.*atan(f/c3)-0.5*c3*log(f.^2+c3^2)),...
    };
N_filters = length(filter_phases_wihtout_chirp);

set(0, 'DefaultLineLineWidth',2)
set(0, 'DefaultAxesColorOrder',[223 0 0]/255)
set(0, 'DefaultAxesFontSize',20)
set(0, 'DefaultFigureRenderer', 'zbuffer');
set(0, 'DefaultFigureWindowStyle', 'normal');
set(0, 'DefaultAxesUnits','pixels')
N_columns = 3;
fig_1_handle = figure('Position',[100 100 550*N_columns 400*(N_filters+1)]);
colors =[255    0    0;
           0  203    0;
           0  173  255]./255;

%% First row plots
fig1_subplot_handles = zeros(N_filters+1, N_columns);
fig1_subplot_handles(1,1) = subplot(N_filters+1,N_columns,1);
plot(t*1e12,abs(E_in_t)/1e3, 'Color', colors(1,:))
xlabel('Time [ps]'), ylabel('Amplitude [kV/m]')
xlim([-320 2200])
ylim([-0.1 0.8])

fig1_subplot_handles(1,2) = subplot(N_filters+1,N_columns,2);
plot(f/1e9,abs(E_in_f)*1e9, 'Color', colors(1,:))
xlabel('Envelope Frequency [GHz]'), ylabel('Magnitude [nV/m/Hz]')
safe_range_of_freq = [-2850 2850].*1e9;
xlim(safe_range_of_freq/1e9)
set(gca,'Xtick',[-2000 0 2000])
ylim([0 11.019])

fig1_subplot_handles(1,3) = subplot(N_filters+1,N_columns,3);
[~,T,F,P] = spectrogram(abs(E_in_f)*1e3, floor(length(E_in_f)*0.0512),...
    floor(length(E_in_f)*0.0500), linspace(0,12e-12,201), 1/df);
design_F = F-Fs/2;
surf(T/1e-12, design_F/1e9, P.', 'edgecolor','none')
colormap hot
view(90,-90)
xlabel 'Freq. of Spectrum [ps]', ylabel 'Envelope Frequency [GHz]'
ylabel(colorbar,'Spectral Density [mV^2/m^2/Hz]','FontSize',15)
axis tight
ylim(safe_range_of_freq/1e9)
set(gca,'Ytick',[-2000 0 2000])
set(gca, 'Xtick',[0 12])
hold on

P_threshold = 0.001;
above_P_threshold = P_threshold.*max(max(P)) < P;
design_effective_bandwidth = zeros(size(F));
for F_index = 1:length(F)
    for T_index = length(T):-1:1
        if above_P_threshold(T_index, F_index)
            design_effective_bandwidth(F_index) = T(T_index);
            break
        end
    end
end
design_delta_design_F_total = 0.5./design_effective_bandwidth;
design_PD_BW = 14.5e9; % Photodetector bandwidth in Hz
design_ADC_BW = 14.5e9; % ADC Nyquist bandwidth in Hz
design_filter_GDD_design_F_DFT = 1./(pi*design_delta_design_F_total.^2);
design_filter_GDD_design_F_PD = 0.35./(2*pi*design_PD_BW*design_delta_design_F_total); 
design_filter_GDD_design_F_ADC = 0.5./(2*pi*design_ADC_BW*design_delta_design_F_total);
design_filter_GDD_design_F = max([design_filter_GDD_design_F_DFT; design_filter_GDD_design_F_PD; design_filter_GDD_design_F_ADC]);
design_dF = mean(diff(design_F));
shifted_design_filter_GD_design_F = cumsum(design_filter_GDD_design_F).*(2*pi*design_dF);
design_filter_GD_design_F = shifted_design_filter_GD_design_F - shifted_design_filter_GD_design_F(value_finder(design_F,0)); % Just to force zero freq. GD to zero
design_filter_GD_f = interp1(design_F, design_filter_GD_design_F, f, 'linear', 'extrap');

chirp_from_filtered_E_in_design_F = zeros(size(F));
for F_index = 1:length(F)
    filter = ((design_F(F_index) - 0.5*design_dF) < f) & (f < (design_F(F_index) + 0.5*design_dF));
    filtered_E_in_t = Fs * ifft(ifftshift(E_in_f.*filter));
    envelope_max_indices = find(abs(filtered_E_in_t) == max(abs(filtered_E_in_t)));
    chirp_from_filtered_E_in_design_F(F_index) = -t(round(mean(envelope_max_indices)));
end
chirp_from_filtered_E_in_f = interp1(design_F, chirp_from_filtered_E_in_design_F, f, 'linear', 'extrap');
chirp_phase_from_filtered_E_in_f = 2*pi*df.*cumsum(chirp_from_filtered_E_in_f);

design_filter_GD_f_chirp_corrected = design_filter_GD_f - chirp_from_filtered_E_in_f;
design_filter_phase_f = 2*pi*df.*cumsum(design_filter_GD_f);
% These are just to show that other methods doesn't work, even Hilbert
% filtering gives instantenous frequency but not the chirp
% normalized_E_in_f = E_in_f./abs(E_in_f);
% normalized_E_in_f_middle = (normalized_E_in_f(2:end) + normalized_E_in_f(1:(end-1)))/2;
% chirp_from_normalized_E_in_f = abs((diff(normalized_E_in_f)/df)./(1j*2*pi*normalized_E_in_f_middle));
% smoothed_chirp_from_normalized_E_in_f = ideal_moving_average(chirp_from_normalized_E_in_f, round(design_df/df));
% chirp_from_angle_E_in_f = diff(angle(E_in_f))/(2*pi*df);
% smoothed_chirp_from_angle_E_in_f = ideal_moving_average(chirp_from_angle_E_in_f, round(design_df/df));
switch 'Design'
    case 'Design'
        filter_phases_wihtout_chirp{2} = @(f) design_filter_phase_f.*f.^0;
        % For testing purposes you can uncomment the following line to show the designed effective bandwidth
        % plot(design_effective_bandwidth/1e-12, design_F/1e9, '-.', 'Color', [255 53 120]/255) %plotting design effective bandwidth
    case 'Do not design'
end

switch 'Unchirp'
    case 'Unchirp'
        filter_phases = cell(1,N_filters);
        for filter_index = 1 : N_filters
            filter_phases{filter_index} = @(f) filter_phases_wihtout_chirp{filter_index}(f) - chirp_phase_from_filtered_E_in_f.*f.^0;
        end
    case 'Do not unchirp'
        filter_phases = filter_phases_wihtout_chirp;
end

fig_2_handle = figure('Position',[100 100 550*2 400*N_filters]);
fig2_subplot_handles = zeros(N_filters, 2);

%% Second and third row plots
for filter_index = 1 : N_filters
    filter_phase_f = filter_phases{filter_index}(f);
    filter_GD_f = diff(filter_phase_f)./(2*pi*df); % Filter group delay [s]
    filter_GD_f_without_chirp = diff(filter_phases_wihtout_chirp{filter_index}(f))./(2*pi*df);
    filter_GDD_f_without_chirp = diff(filter_GD_f_without_chirp)./(2*pi*df); % Filter group delay dispersion [s^2]
    filter_f = exp(1j*filter_phase_f); % Complex filter: down-converted filter
    E_out_f = filter_f.*E_in_f; % Complex envelope of output electric field in frequency [V/m/Hz]
    E_out_t = Fs * ifft(ifftshift(E_out_f)); % Complex envelope of output electric field in time [V/m]
    E_out_envelope_t = abs(E_out_t); % Envelope of output electric field in time [V/m]
    E_out_envelope_f = 1/Fs * fftshift(fft(E_out_envelope_t)); % Envelope of output electric field in frequency [V/m/Hz]
    conversion_factor = c*3.5*8.854e-12/2*(100e-6)^2*1*100; % Conversion equation: (c*n*e0/2*(abs(E_t))^2)*A*r*G [A: area, r: responsivity, G: gain]
    i_out_envelope_t = conversion_factor*(abs(E_out_t)).^2; % Envelope of output photocurrent in time [A]
    i_out_envelope_f = 1/Fs * fftshift(fft(i_out_envelope_t)); % Envelope of output photocurrent in frequency [A/Hz]
    switch filter_index
        case 1
            PD_BW = 14.5e9; % Photodetector bandwidth in Hz
            ADC_BW = 14.5e9; % ADC Nyquist bandwidth in Hz
        case 2
            PD_BW = 14.5e9; % Photodetector bandwidth in Hz
            ADC_BW = 14.5e9; % ADC Nyquist bandwidth in Hz
    end
    Aquisition_BW = min(PD_BW, ADC_BW); % Minimum of all electrical bandwidth limitations
    delta_f_DFT = sqrt(1./(pi*abs(filter_GDD_f_without_chirp))); % Spectral resolution limit imposed by Dispersive Fourier Transform
    delta_f_PD = 0.35./(2*pi*PD_BW*abs(filter_GDD_f_without_chirp)); % Spectral resolution limit imposed by Photodetector bandwidth
    delta_f_ADC = 0.5./(2*pi*ADC_BW*abs(filter_GDD_f_without_chirp)); % Spectral resolution limit imposed by ADC Nyquist bandwidth
    delta_f_total = max([delta_f_DFT; delta_f_PD; delta_f_ADC]); % Overall spectral resolution limit
    Largest_freq_of_spectrum = 0.5./delta_f_total; % Effective bandwidth of spectrum's spectrum
    
    figure(fig_2_handle)
    fig2_subplot_handles(filter_index, 1) = subplot(N_filters,2,2*(filter_index-1)+1);
    plot(f_GDD./1e9, delta_f_DFT/1e9, '-', 'Color', [0.85 0.325 0.098])
    hold on
    plot(f_GDD./1e9, delta_f_PD/1e9, '-', 'Color', [0.929 0.694 0.125])
    plot(f_GDD./1e9, delta_f_ADC/1e9, '-', 'Color', [0.494 0.184 0.556])
    plot(f_GDD./1e9, delta_f_total/1e9, '--', 'Color', colors(filter_index+1,:))
    xlabel('Envelope Frequency [GHz]'), ylabel('Resolution [GHz]')
    xlim(safe_range_of_freq/1e9)
%     ylim([-10 600])
    legend_resolution_handle = legend('DFT', 'PD', 'ADC', 'Total');
    set(legend_resolution_handle, 'Location', 'NorthEast')
    set(legend_resolution_handle, 'Box', 'on')
    set(legend_resolution_handle, 'FontSize', 15)
    
    fig2_subplot_handles(filter_index, 2) = subplot(N_filters,2,2*(filter_index-1)+2);
    within_safe_range_of_freq = (safe_range_of_freq(1)<f) & (f<safe_range_of_freq(2));
    everyother = 1; % Determines the distance between consecuitive points that the errorbar is plotted for
    magnification_factor = 10; % Magnification factor to visualize the spectral resolution better
    within_safe_range_of_freq_everyother = within_safe_range_of_freq & (mod(1:N,everyother) == 0);
    errorbar_handle = herrorbar((f_GD(within_safe_range_of_freq_everyother(1:end-1))+f_GD(within_safe_range_of_freq_everyother(2:end)))/2/1e9,...
        (filter_GD_f(within_safe_range_of_freq_everyother(1:end-1))+filter_GD_f(within_safe_range_of_freq_everyother(2:end)))/2*1e9,...
        (magnification_factor*delta_f_total(within_safe_range_of_freq_everyother(2:end-1))/2)/1e9);
    set(errorbar_handle(1), 'Color', colors(filter_index+1,:))
    set(errorbar_handle(2), 'Color', colors(filter_index+1,:).^0.2)
    xlabel('Envelope Frequency [GHz]'), ylabel('Group Delay [ns]')
    xlim(safe_range_of_freq/1e9)
    ylim([-3 3])
    set(gca, 'XTick', [-2000 0 2000])
    
    figure(fig_1_handle)
    fig1_subplot_handles(filter_index+1,1) = subplot(N_filters+1,N_columns,N_columns*(filter_index)+1);
    plot(f_GD/1e9, filter_GD_f_without_chirp/1e-9, 'Color', [255 53 120]/255)
    hold on
    plot(f/1e9, chirp_from_filtered_E_in_f/1e-9, 'Color', [255 208 0]/255)
    plot(f_GD/1e9, filter_GD_f/1e-9, 'Color', colors(filter_index+1,:))
    xlabel('Envelope Frequency [GHz]'), ylabel('Group Delay [ns]')
    xlim(safe_range_of_freq/1e9)
    ylim([-3 3])
    set(gca, 'XTick', [-2000 0 2000])
    legend_design_handle = legend('Design', 'Chirp', 'Total');
    set(legend_design_handle, 'Location', 'SouthEast')
    set(legend_design_handle, 'Box', 'on')
    set(legend_design_handle, 'FontSize', 15)
    
    fig1_subplot_handles(filter_index+1,2) = subplot(N_filters+1,N_columns,N_columns*(filter_index)+2);
    plot(t*1e9, abs(E_out_envelope_t), 'Color', colors(filter_index+1,:))
    xlabel('Time [ns]'), ylabel('Amplitude [V/m]')
    xlim([-1.71 1.71])
    ylim([0 450])
    
    fig1_subplot_handles(filter_index+1,3) = subplot(N_filters+1,N_columns,N_columns*(filter_index)+3);
    [S,F,T,P] = spectrogram(E_out_envelope_t*1e3, floor(length(E_out_envelope_t)*0.0512),...
        floor(length(E_out_envelope_t)*0.0500), linspace(0,20e9,201), Fs);
    surf(F/1e9,(T-max_t)*1e9,P.','edgecolor','none')
    colormap hot
    view(90,-90)
    xlabel 'Frequency [GHz]', ylabel 'Time [ns]'
    ylabel(colorbar,'Spectral Density [mV^2/m^2/Hz]','FontSize',15)
    axis tight
    ylim([-1.71 1.71])
    set(gca, 'Xtick',[0 20])
    hold on
    input_spectrum_show_range = 0.99*get(fig1_subplot_handles(1,2),'XLim');
    input_spectrum_GD_starts = filter_GD_f_without_chirp(value_finder(f_GD/1e9, input_spectrum_show_range(1)));
    input_spectrum_GD_ends = filter_GD_f_without_chirp(value_finder(f_GD/1e9, input_spectrum_show_range(2)));
    plot([0 Aquisition_BW Aquisition_BW 0]/1e9, [-input_spectrum_GD_ends(end) -input_spectrum_GD_ends(end) ...
        -input_spectrum_GD_starts(1) -input_spectrum_GD_starts(1)]/1e-9,'-.', 'Color', colors(filter_index+1,:))
    
    axes(fig1_subplot_handles(1,3)) %#ok<LAXES> % This plots on the first row spectrogram
    input_spectrum_GDD_start_indices = value_finder(f_GDD/1e9, input_spectrum_show_range(1));
    input_spectrum_GDD_end_indices = value_finder(f_GDD/1e9, input_spectrum_show_range(2));
    rectangle_indices = [input_spectrum_GDD_start_indices(1),...
        input_spectrum_GDD_start_indices(1):input_spectrum_GDD_end_indices(end),...
        input_spectrum_GDD_end_indices(end)];
    plot([0, Largest_freq_of_spectrum(input_spectrum_GDD_start_indices(1):input_spectrum_GDD_end_indices(end)), 0]/1e-12, f_GDD(rectangle_indices)/1e9,'-.', 'Color', colors(filter_index+1,:))
end
desired_Position = get(fig1_subplot_handles(1,1), 'Position');
for row = 1 : (N_filters+1)
    current_Position = get(fig1_subplot_handles(row,3), 'Position');
    set(fig1_subplot_handles(row,3), 'Position',[current_Position(1:2), desired_Position(3), current_Position(4)]);
end
% Change the colormap
load('modified_hot_colormap.mat')
colormap(fig_1_handle, modified_hot_colormap)

%% Panel letters and arrows
figure(fig_1_handle)
for row = 1 : (N_filters+1)
    for column = 1 : N_columns
        current_Position = get(fig1_subplot_handles(row, column), 'Position');
        annotation('textbox', 'Units','pixels', 'Position',[current_Position(1)-110, current_Position(2)+115, 200, 200],...
            'String', ['(', char(96 + N_columns*(row-1) + column), ')'], 'FontSize',20, 'FontWeight','bold', 'FontName','calibri', 'LineStyle','none');
    end
end

figure(fig_2_handle)
for row = 1 : N_filters
    for column = 1 : 2
        current_Position = get(fig2_subplot_handles(row, column), 'Position');
        annotation('textbox', 'Units','pixels', 'Position',[current_Position(1)-110, current_Position(2)+115, 200, 200],...
            'String', ['(', char(96 + 2*(row-1) + column), ')'], 'FontSize',20, 'FontWeight','bold', 'FontName','calibri', 'LineStyle','none');
    end
    switch row
        case 1
            arrow_tips_distance = 84; % in pixels 
        case 2
            arrow_tips_distance = 84; % in pixels
    end
    annotation('arrow', 'Units','pixels', 'Position',[current_Position(1)+(current_Position(3)/2)-arrow_tips_distance, current_Position(2)+(current_Position(4)/2), 60, 0]);
    annotation('arrow', 'Units','pixels', 'Position',[current_Position(1)+(current_Position(3)/2)+arrow_tips_distance, current_Position(2)+(current_Position(4)/2), -60, 0]);
    annotation('textbox', 'Units','pixels', 'Position',[current_Position(1)+(current_Position(3)/2)-arrow_tips_distance-20, current_Position(2)+(current_Position(4)/2)-20, 100, 100],...
        'String', [num2str(magnification_factor), ' \times Resolution'], 'FontSize',15, 'FontWeight','bold', 'FontName','calibri', 'LineStyle','none', 'HorizontalAlign','right');
end
