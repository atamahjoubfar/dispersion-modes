lambda = 1.5e-6; % Carrier wavelength [m]
c = 299792458; % Speed of light [m/s]
f0 = c/lambda; % Carrier frequency [Hz]

Fs = 1e13; % Sampling rate [Hz]
dt = 1/Fs; % Time simulation resolution [s]
max_t = 1e-8; % Maximum simulation time [s]
t = -max_t : dt : max_t; % Time samples [s]
N = length(t); % Number of samples
df = 0.5/max_t; % Frequency simulation resolution [Hz]
f = -0.5/dt : df : 0.5/dt; % Frequency samples [Hz]

switch 'set temporal shape'
    case 'set temporal shape'
        E_in_t = 1e3*(exp(-((t-4e-12)/1e-12).^2)+0.7*exp(-((t+4e-12)/1e-12).^2)); % Complex envelope of input electric field in time [V/m]
        E_in_f = 1/Fs * fftshift(fft(E_in_t)); % Complex envelope of input electric field in frequency [V/m/Hz]
    case 'set spectral shape'
        E_in_f = exp(-((f-4e10)/1e10).^2)+exp(-((f+4e10)/1e10).^2); %.*exp(1j*pi*rand(1,N));
        E_in_t = Fs * ifft(ifftshift(E_in_f));
end
E_in_envelope_t = abs(E_in_t); % Envelope of input electric field in time [V/m]
E_in_envelope_f = 1/Fs * fftshift(fft(E_in_envelope_t)); % Envelope of input electric field in frequency [V/m/Hz]

c0 = 2*pi*300*(1e-12).^0;
c1 = 2*pi*300*(1e-12).^1;
c2 = 2*pi*300*(1e-12).^2;
c3 = 2*pi*300*(1e-12).^3;
c4 = 2*pi*300*(1e-12).^4;
c5 = 2*pi*300*(1e-12).^5;
c6 = 2*pi*300*(1e-12).^6;
filter_phases = { @(f) c0.*f.^0,...
    @(f) c1.*f.^1,...
    @(f) c2.*f.^2,...
    @(f) c3.*f.^3,...
    @(f) c4.*f.^4,...
    @(f) c5.*f.^5,...
    @(f) c6.*f.^6,...
    };
N_filters = length(filter_phases);

set(0, 'DefaultLineLineWidth',1)
set(0, 'DefaultAxesFontSize',12)
set(0, 'DefaultFigureRenderer', 'zbuffer');
set(0, 'DefaultFigureWindowStyle', 'normal');
figure('Position', [1, 1, 2*1169, 1513])
columns_to_show = {'Phase', 'GD', 'GDD', '|E_t|', '|F|E_t||', 'i_t', 'i_f'}; % {'Phase', 'GD', 'GDD', 'E_t', 'E_f', 'E_spectrogram', 'i_t', 'i_f', 'i_spectrogram', 'Resolution', 'Effective_BW'}
N_columns = length(columns_to_show);
handles = zeros(N_filters, N_columns);

for filter_index = 1 : N_filters
    filter_phase_f = filter_phases{filter_index}(f);
    filter_GD_f = diff(filter_phase_f)./(2*pi*df); % Filter group delay [s]
    f_GD = (f(2:end)+f(1:end-1))/2; % Frequencies of GD samples
    filter_GDD_f = diff(filter_GD_f)./(2*pi*df); % Filter group delay dispersion [s^2]
    f_GDD = f(2:end-1); % Frequencies of GDD samples
    filter_f = exp(1j*filter_phase_f); % Complex filter: down-converted filter
    E_out_f = filter_f.*E_in_f; % Complex envelope of output electric field in frequency [V/m/Hz]
    E_out_t = Fs * ifft(ifftshift(E_out_f)); % Complex envelope of output electric field in time [V/m]
    E_out_envelope_t = abs(E_out_t); % Envelope of output electric field in time [V/m]
    E_out_envelope_f = 1/Fs * fftshift(fft(E_out_envelope_t)); % Envelope of output electric field in frequency [V/m/Hz]
    conversion_factor = c*3.5*8.854e-12/2*(100e-6)^2*1*100; % Conversion equation: (c*n*e0/2*(abs(E_t))^2)*A*r*G [A: area, r: responsivity, G: gain]
    i_out_envelope_t = conversion_factor*(abs(E_out_t)).^2; % Envelope of output photocurrent in time [A]
    i_out_envelope_f = 1/Fs * fftshift(fft(i_out_envelope_t)); % Envelope of output photocurrent in frequency [A/Hz]
    PD_BW = 80e9; % Photodetector bandwidth in Hz
    ADC_BW = 120e9; % ADC Nyquist bandwidth in Hz
    Aquisition_BW = min(PD_BW, ADC_BW); % Minimum of all electrical bandwidth limitations
    delta_f_DFT = sqrt(abs(1./(pi*filter_GDD_f))); % Spectral resolution limit imposed by Dispersive Fourier Transform
	delta_f_PD = 0.35./(2*pi*PD_BW*filter_GDD_f); % Spectral resolution limit imposed by Photodetector bandwidth
	delta_f_ADC = 0.5./(2*pi*ADC_BW*filter_GDD_f); % Spectral resolution limit imposed by ADC Nyquist bandwidth
	delta_f_total = max([delta_f_DFT; delta_f_PD; delta_f_ADC]); % Overall spectral resolution limit
    Largest_freq_of_spectrum = 0.5./delta_f_total; % Effective bandwidth of spectrum's spectrum

    for column_index = 1 : N_columns
        handles(filter_index,column_index) = subplot(N_filters,N_columns,N_columns*(filter_index-1)+column_index);
        switch columns_to_show{column_index}
            case 'Phase'
                plot(f/1e9, filter_phase_f, 'Color',[255 0 102]/255, 'LineWidth',2)
                xlabel('Frequency [GHz]'), ylabel('\phi(f) [Radians]')
                
            case 'GD'
                plot(f_GD/1e9, filter_GD_f.*1e9, 'Color',[0 204 0]/255, 'LineWidth',2)
                xlabel('Frequency [GHz]'), ylabel('GD(f) [ns]')
                
            case 'GDD'
                plot(f_GDD/1e9, filter_GDD_f.*(1e9)^2, 'Color',[31 111 255]/255, 'LineWidth',2)
                xlabel('Frequency [GHz]'), ylabel('GDD(f) [ns^2]')
                
            case '|E_t|'
                plot(t*1e9, E_out_envelope_t, 'Color',[255 128 0]/255, 'LineWidth',0.5)
                axisp(99.9, 'tight')
                xlabel('Time [ns]'), ylabel('|E(t)| [V/m]')
                
            case '|F|E_t||'
                plot(f/1e9, abs(E_out_envelope_f), 'Color',[255 128 0]/255, 'LineWidth',0.5)
                axisp(99.9, 'tight')
                xlabel('Frequency [GHz]'), ylabel('|F\{|E(t)|\}| [V/m/Hz]')
                
            case '|E_f|'
                plot(f/1e9, abs(E_out_f), 'Color',[255 128 0]/255, 'LineWidth',0.5)
                axisp(99.9, 'tight')
                xlabel('Frequency [GHz]'), ylabel('|E(f)| [V/m/Hz]')
                
            case 'E_spectrogram'
                [S,F,T,P] = spectrogram(E_out_envelope_t,256,250,256,Fs);
                surf(T*1e9,F/1e9,P,'edgecolor','none')
                axis tight, view(90,-90)
                xlabel 'Time [ns]', ylabel 'Frequency [GHz]'
                
            case 'i_t'
                plot(t*1e9, i_out_envelope_t, 'Color',[122 16 228]/255, 'LineWidth',0.5)
                axisp(99.9, 'tight')
                xlabel('Time [ns]'), ylabel('I(t) [A]')
                
            case 'i_f'
                plot(f/1e9, abs(i_out_envelope_f), 'Color',[122 16 228]/255, 'LineWidth',0.5)
                axisp(99.9, 'tight')
                xlabel('Frequency [GHz]'), ylabel('I(f) [A/Hz]')
                
            case 'i_spectrogram'
                [S,F,T,P] = spectrogram(i_out_envelope_t,256,250,256,Fs);
                surf(T*1e9,F/1e9,P,'edgecolor','none')
                axis tight, view(90,-90)
                xlabel 'Time [ns]', ylabel 'Frequency [GHz]'

            case 'Effective_BW'
                plot(f_GDD/1e9, Largest_freq_of_spectrum/1e-9, 'Color',[122 122 122]/255, 'LineWidth',0.5)
                xlabel('Frequency [GHz]'), ylabel('BW_{spectrum}(f) [ns]')

            case 'Resolution'
           		plot(f_GDD./1e9, delta_f_total/1e9, 'Color',[122 122 122]/255, 'LineWidth',0.5)
				xlabel('Frequency [GHz]'), ylabel('Resolution [GHz]')
        end
    end
end

set(handles(1,2),'Ylim',[-1 1])
set(handles(2,2),'Ylim',[-1 1])
