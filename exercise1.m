%% analyzing a sound file - assignment 1
clc, close all; clear all;
% read sound file
[audio, sample_rate] = audioread("tone_in_noise.wav");

% plot sound file
signal_length = numel(audio);
timeseries = linspace(0, signal_length/sample_rate, signal_length);

figure('Position', [100, 100, 800, 400]);
plot(timeseries, audio)
xlabel("Time (s)")
ylabel("Signal amplitude (Pa)")
xlim([0 max(timeseries)])
set(gca, 'FontSize', 14);
grid on;

% calculate and plot single-sided amplitude spectrum
[fft_positive, frequencies] = single_sided_fft(audio, sample_rate);

figure('Position', [100, 100, 800, 400]);
loglog(frequencies, fft_positive, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (Pa)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

% -> single tone: frequency = 3000Hz and amplitude = 0.00459732

%% assignment 2
% cut of heading and trailing zeros in the signal
lower_cutoff = sample_rate * 0.7;
upper_cutoff = sample_rate * 1.7;
audio_cut = audio(lower_cutoff:upper_cutoff);

% plot cut audio
% signal_length_cut = numel(audio_cut);
% timeseries_cut = linspace(0, signal_length_cut/sample_rate, signal_length_cut);
% plot(timeseries_cut, audio_cut)
% xlabel("Time (seconds)")
% ylabel("Signal amplitude (Pa)")
% set(gca, 'FontSize', 14);
% grid on;

% calculate and plot single-sided amplitude spectrum
[fft_positive, frequencies] = single_sided_fft(audio_cut, sample_rate);

figure('Position', [100, 100, 800, 400]);
loglog(frequencies, fft_positive, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (Pa)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

% -> single tone: frequency = 2999,94Hz and amplitude = 0.00996987

%% generation of a signal - assignment 3
clc, close all, clear all;

% create timeseries
sample_rate = 10240;
samples = 1024;

% generate signal
times = (0:samples)'/sample_rate;
signal1 = 0.5 + 0.1 * sin(2 * pi * 100 * times);

% calculate and plot amplitude spectrum of signal1
[fft_positive, frequencies] = single_sided_fft(signal1, sample_rate);

figure('Position', [100, 100, 400, 400]);
loglog(frequencies, fft_positive, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (V)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

% generate signal but with modified timepoints (signal1_mod)
times_mod = (0:samples-1)'/sample_rate;
signal1_mod = 0.5 + 0.1 * sin(2 * pi * 100 * times_mod);

% calculate and plot amplitude spectrum of signal1_mod
[fft_positive, frequencies] = single_sided_fft(signal1_mod, sample_rate);

figure('Position', [100, 100, 400, 400]);
loglog(frequencies, fft_positive, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (V)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

%% assignment 4

% create timeseries
sample_rate = 10240;
samples = 1024;
times_mod = (0:samples-1)'/sample_rate;
% generate signal 
signal2 = 0.1 * sin(2 * pi * 105 * times_mod);

% calculate rms
signal2_rms = rms(signal2)

% -> signal2_rms = 0.0707

% calculate and plot amplitude spectrum of signal2
[fft_positive, frequencies] = single_sided_fft(signal2, sample_rate);

figure('Position', [100, 100, 400, 400]);
loglog(frequencies, fft_positive, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (V)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

% -> tone amplitude = 0.065149 and frequency=99.9024

% calculate signal amplitude maximum error = signal2_rms - tone amplitude
% -> max_error = 0.005551

%% assignment 5
% create flat top window
window_length = samples;
flattop_window = flattopwin(window_length);

% apply flat top window
signal2_windowed = signal2 .* flattop_window * samples/sum(flattop_window);
rms_signal2windowed = rms(signal2_windowed)

% plot amplitude spectrum of signal2*window
[fft_positive, frequencies] = single_sided_fft(signal2_windowed, sample_rate);

figure('Position', [100, 100, 500, 400]);
loglog(frequencies, fft_positive, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (V)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

%% noise and power spectra
%% assignment 6
clc, close all, clear all;
% generate signal
sample_rate3 = 1024;
duration = 1;
signal_length3 = sample_rate3 * duration;
times3 = (0:signal_length3-1)/sample_rate3;

signal3 = 0.5 + 0.1 * sin(2 * pi * 100 * times3);

% generating two more signals
% define signal to noise ratio in dB
snr = 0;

signal4 = awgn(signal3, snr);
signal5 = awgn(signal3, snr);

% plot the waveform of the signals
figure;
subplot(3,1,1);
plot(times3, signal3);
title('Signal 3');
xlabel('Time (s)');
ylabel('Signal amplitude (V)');
xlim([0, max(times3)]);
set(gca, 'FontSize', 14);
grid on;

subplot(3,1,2);
plot(times3, signal4);
title('Signal 4');
xlabel('Time (s)');
ylabel('Signal amplitude (V)');
xlim([0, max(times3)]);
set(gca, 'FontSize', 14);
grid on;

subplot(3,1,3);
plot(times3, signal5);
title('Signal 5');
xlabel('Time (s)');
ylabel('Signal amplitude (V)');
xlim([0, max(times3)]);
set(gca, 'FontSize', 14);
grid on;

% calculate power spectrum of signal4 and 5
[fft_positive, frequencies] = single_sided_fft(signal4, sample_rate3);
fft_signal4 = fft_positive;
[fft_positive, frequencies] = single_sided_fft(signal5, sample_rate3);
fft_signal5 = fft_positive;

% plot amplitude spectrum of signal4 and 5
figure('Position', [100, 100, 500, 400]);
loglog(frequencies, fft_signal4, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (V)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

figure('Position', [100, 100, 500, 400]);
loglog(frequencies, fft_signal5, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Signal Amplitude (V)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

%calculate power spectra
powerspectra4 = abs(fft_signal4) .^ 2;
powerspectra5 = abs(fft_signal5) .^ 2;

% convert to dB
powerspectra4_dBm = 10 * log10(powerspectra4/(600e-3));
powerspectra5_dBm = 10 * log10(powerspectra5/(600e-3));

% plot powerspectra of signal4
figure('Position', [100, 100, 500, 400]);
plot(frequencies, powerspectra4, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Power Spectra (dBm)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

% plot powerspectra of signal5
figure('Position', [100, 100, 500, 400]);
plot(frequencies, powerspectra5_dBm, "LineWidth", 1) 
xlabel("Frequency (Hz)")
ylabel("Power Spectra (dBm)")
xlim([min(frequencies) max(frequencies)])
set(gca, 'FontSize', 14);
grid on;

%% assignment 5 continuation (execute section before to avoid errors)
% calculate power from time domain
% calculate rms
signal4_rms = rms(signal4)
signal5_rms = rms(signal5)
% signal4_rms2 = sqrt(1 / signal_length3 * sum(signal4.^2)); % same values
% calculate power from rms
power_rms4 = signal4_rms ^ 2
power_rms5 = signal5_rms ^ 2
% note:
    % power_rms has the unit Volts squared
    % do not divide by the resistance because then the unit is Watts
    % but the unit of the power_ftt (see below) is only Volts squared
    % if you divide by the resistance the order of magnitude of both values
    % will not be correct

% caluclate power from frequency domain
% sum power spectra
power_fft4 = sum(powerspectra4)
power_fft5 = sum(powerspectra5)

%% sound file revisited - assignment 7
% no programming

%% play sound - assignment 8
clc, close all; clear all;

% define variables
sample_rate = 48e3;
duration = 0.5;
signal_length = sample_rate * duration;
times = (0:signal_length-1)'/sample_rate;
% generate signal
signal6 = 0.1 * cos(2 * pi * 230 * times);

% generate ramp filter with appropriate signal_length
ramp_time = 100e-3;
ramp_length = sample_rate * ramp_time;
ramp_filter = create_ramp_filter(ramp_length, signal_length);
% plot(ramp_filter)

% apply ramp filter to signal
signal6_rampfiltered = signal6 .* ramp_filter;

% save rampfiltered signal6
audiowrite("signal6_rampfiltered.wav", signal6_rampfiltered, sample_rate)

figure('Position', [100, 100, 500, 400]);
plot(times, signal6_rampfiltered)
hold on
plot(times, 0.1*ramp_filter, 'LineWidth', 3, 'Color', 'r')
legend('Signal E', "Ramp Filter", 'Location', 'southeast')
xlabel("Time (s)")
ylabel("Signal amplitude(V)")
xlim([0 max(times)])
set(gca, 'FontSize', 14);
grid on;



