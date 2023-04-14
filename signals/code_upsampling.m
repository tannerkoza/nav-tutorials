%% Code Upsampling Tutorial
%This script outlines upsampling codes (eg. PRNs) for linearly changing Doppler
%frequencies.

clear
clc
close all

%% User-Defined Varaibles

% time & sampling
sig_duration = 0.005; % signal duration [s]
fsamp = 40000; % sampling frequency [Hz]
num_code_periods = 10; % integration duration [s]

% Doppler
max_doppler = 50; % [Hz]

% baseband code
fcode = 10; % code frequency [Hz]
fchip = 100; % chip frequency [Hz or chips/s]

%% Initialization

% Doppler
doppler_profile = linspace(0, max_doppler, num_code_periods); % [Hz or chips/s]

% baseband code
samp_per_chip = fsamp / fchip; % samples per chip
chip_per_code = fchip / fcode; % chips per code period
rem_phase = 0;

code = 2 * randi([0, 1], 1, chip_per_code) - 1; % NRZ code

% plotting
min_samp_idx = 0;

%% Upsampling

figure('Name','Upsampled Code with Doppler Shift')
hold on

for per = 1:num_code_periods
    % upsampling
    fchip = fchip + doppler_profile(per);
    [upsamp_code, rem_phase] = upsample(code, rem_phase, fsamp, fchip);

    % plotting
    max_samp_idx = min_samp_idx + length(upsamp_code)-1; % subtract 1 for zero-numbering
    samp_idx = min_samp_idx:max_samp_idx;
    min_samp_idx = max_samp_idx + 1;
    plot(samp_idx,upsamp_code,'b')
    xline(max_samp_idx,'r','LineWidth',2)
end

title('Upsampling Example with Code Doppler')
xlabel('# of Samples')
legend('Upsampled Code','Code Period End','Location','best')

function [upsamp_code, new_rem_phase] = upsample(code, rem_phase, fsamp, fchip)
        % initialization
        code_length = length(code);
        samp_per_chip = fsamp/fchip;
        chip_per_samp = 1/samp_per_chip;
        samp_per_code_period = ceil((code_length-rem_phase)/chip_per_samp);
        appended_code = [code(end) code code(1)];

        % upsampling
        code_subchip_idx = rem_phase:chip_per_samp:(samp_per_code_period-1)*chip_per_samp + rem_phase; % [fractional chips]
        code_chip_idx = ceil(code_subchip_idx) + 1; % add 1 for one-indexing [whole chips]
        upsamp_code = appended_code(code_chip_idx);

        % remainder code phase calculation
        new_rem_phase = code_subchip_idx(samp_per_code_period) + chip_per_samp - code_length; % [fractional chips]
end
