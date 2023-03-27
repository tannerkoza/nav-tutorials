%% Phase Shift Keying Tutorial
%This script includes basic BPSK and QPSK modulation and demodulation
%schemes. The acqusition of these signals using the modulated code is also
%conducted.

clear
clc
close all

%% User-Defined Variables

% time & sampling
sig_duration = 1.1; % signal duration [s]
fsamp = 10000; % sampling frequency [Hz]

% carrier signal
fcarr = 10; % carrier frequency [Hz]

% baseband code
fcode = 10; % code frequency [Hz]
fchip = 1000; % chip frequency [Hz]
code_data_shift = 45; % code shift [chips]
code_pilot_shift = 20;
code_per_duration = sig_duration/(1/fcode);

%% Initialization

% time & sampling
t = 0:(1 / fsamp):(sig_duration - 1 / fsamp); % time vector [s]

% carrier signal
carr = exp(1i*2*pi*fcarr*t); % carrier signal
carr_rep = exp(-1i*2*pi*fcarr*t); % carrier signal replica (negative frequency)

% baseband code
samp_per_chip = fsamp / fchip; % samples per chip
chip_per_code = fchip / fcode; % chips per code period

code_data = 2 * randi([0, 1], 1, chip_per_code) - 1; % NRZ data channel code
upsamp_code_data = repvec(repelem(code_data, samp_per_chip), code_per_duration);
shift_code_data = circshift(upsamp_code_data, code_data_shift*samp_per_chip); 

code_pilot = 2 * randi([0, 1], 1, chip_per_code) - 1; % NRZ pilot channel code
upsamp_code_pilot = repvec(repelem(code_pilot, samp_per_chip), code_per_duration);
shift_code_pilot = circshift(upsamp_code_pilot, code_pilot_shift*samp_per_chip); 

%% BPSK Modulation
%To simulate a real signal, only use the real part of the carrier.

bpsk = carr .* shift_code_data; % IQ modulation

figure('Name','BPSK Modulation')
tiledlayout(2, 1)
nexttile
plot(t, real(bpsk))
title('In-Phase BPSK Signal')

nexttile
plot(t, imag(bpsk))
title('Quadra-Phase BPSK Signal')
xlabel('time [s]')

%% BPSK Demodulation
%To demodulate real signal, low pass filter in-phase of baseband or
%subtract replicated signal with code (not practical).

baseband_bpsk = bpsk .* carr_rep;
[bpsk_correlation, bpsk_sample_lag] = acquire(baseband_bpsk, upsamp_code_data);
plot_acquisition(bpsk_correlation)

figure('Name','BPSK Demodulation')
tiledlayout(2, 1)
nexttile
hold on
plot(t, real(baseband_bpsk),'DisplayName','Wiped Signal')
plot(t, 0.5*shift_code_data,'DisplayName','Replica Data Code')
title('In-Phase Wiped BPSK Signal')
legend

nexttile
plot(t, imag(baseband_bpsk))
title('Quadra-Phase Wiped BPSK Signal')
xlabel('time [s]')

%% QPSK Modulation

qpsk_code = shift_code_data + shift_code_pilot*1i;
qpsk = carr .* qpsk_code;

figure('Name','QPSK Modulation')
tiledlayout(2, 1)
nexttile
plot(t, real(qpsk))
title('In-Phase QPSK Signal')

nexttile
plot(t, imag(qpsk))
title('Quadra-Phase QPSK Signal')
xlabel('time [s]')

%% QPSK Demodulation
 
baseband_qpsk = qpsk.*carr_rep;

[qpsk_correlation, qpsk_sample_lag] = acquire(baseband_qpsk, upsamp_code_pilot);
plot_acquisition(qpsk_correlation)

figure('Name','QPSK Demodulation')
tiledlayout(2, 1)
nexttile
hold on
plot(t, real(baseband_qpsk),'DisplayName','Wiped Signal')
plot(t, 0.5*shift_code_data,'DisplayName','Replica Data Code')
title('In-Phase Wiped QPSK Signal')
legend

nexttile
hold on
plot(t, imag(baseband_qpsk),'DisplayName','Wiped Signal')
plot(t, 0.5*shift_code_pilot,'DisplayName','Replica Pilot Code')
title('Quadra-Phase Wiped QPSK Signal')
xlabel('time [s]')
legend

%% Acquisition

function [correlation, sample_lag] = acquire(baseband_sig, replica)
    % FFT
    baseband_fft = fft(baseband_sig);
    replica_fft = fft(replica);

    % Correlate
    correlation_fft = baseband_fft .* conj(replica_fft);
    correlation = abs(ifft(correlation_fft)).^2;
    
    [~, correlation_idx] = max(correlation);
    sample_lag = correlation_idx - 1;
end

function plot_acquisition(correlation)
    num_samp = length(correlation);
    samp_idx = 0:num_samp-1;

    figure('Name','PSK Acquisition')
    plot(samp_idx, correlation)
    xlabel('Sample Index')
    ylabel('Correlation Magnitude')
    axis padded
    title('PSK Acquistion Correlation')
end

function out_vec = repvec(vec, mult)
    out_vec = [repmat(vec,1,fix(mult)) vec(1:round(numel(vec)*mod(mult,1)))];
end