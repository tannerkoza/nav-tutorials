%% Carrier Tracking Tutorial

clear
clc
close all

%% User-Defined Variables

% time & sampling
duration = 5; % [s]
int_period = 0.001; % [s]

% signal
fcarrier = 1000; % [Hz]
fcarrierrate = 50; %[Hz/s]
fsamp = 1e6; % [Hz]
fbasis = 960; % [Hz]
cn0 = 60; % [dBHz]

% filter settings
fll_noise_bw = 10; % [Hz]
pll_noise_bw = 15; % [Hz]

%% Initialization

% time & sampling
samp_period = 1 / fsamp;
time = 0:int_period:(duration - int_period);
num_int_periods = ceil(duration/int_period);

% signal
carr = navtools.dsp.complex_carrier(fcarrier, fcarrierrate, fsamp, duration, cn0);
start = 1;
stop = 1 + fsamp * int_period - 1;

% locked loops
fll_ip_last = 0;
fll_qp_last = 0;
fll_w = 0;
fll_x = 0;
pll_w = 0;
pll_x = 0;
fll_rem_phase = 0;
fll_replica_freq = fbasis;
pll_rem_phase = 0;
pll_replica_freq = fbasis;

for period = 1:num_int_periods
    carr_in_period = carr(start:stop);
    [fll_replica, fll_rem_phase] = create_carr_replica(fll_replica_freq, fsamp, fll_rem_phase, length(carr_in_period));
    [pll_replica, pll_rem_phase] = create_carr_replica(pll_replica_freq, fsamp, pll_rem_phase, length(carr_in_period));

    fll_ip = real(sum(carr_in_period.*fll_replica));
    fll_qp = imag(sum(carr_in_period.*fll_replica));
    pll_ip = real(sum(carr_in_period.*pll_replica));
    pll_qp = imag(sum(carr_in_period.*pll_replica));

    fll_error = freq_disc(fll_ip, fll_ip_last, fll_qp, fll_qp_last, int_period);
    pll_error = phase_disc(pll_ip, pll_qp);

    [fll_freq_offset, fll_w, fll_x] = freq_lock_loop(fll_error, fll_noise_bw, fll_w, fll_x, int_period);
    [pll_freq_offset, pll_w, pll_x] = phase_lock_loop(pll_error, pll_noise_bw, pll_w, pll_x, int_period);

    fll_replica_freq = fbasis + fll_freq_offset;
    pll_replica_freq = fbasis + pll_freq_offset;

    start = start + fsamp * int_period;
    stop = stop + fsamp * int_period;

    fll_ip_last = fll_ip;
    fll_qp_last = fll_qp;

    fll_ip_log(period) = fll_ip;
    fll_qp_log(period) = fll_qp;
    fll_error_log(period) = fll_error;
    fll_replica_freq_log(period) = fll_replica_freq;
    fll_freq_offset_log(period) = fll_freq_offset;

    pll_ip_log(period) = pll_ip;
    pll_qp_log(period) = pll_qp;
    pll_error_log(period) = pll_error;
    pll_replica_freq_log(period) = pll_replica_freq;
    pll_freq_offset_log(period) = pll_freq_offset;

end

figure('Name','FLL Error')
plot(time, fll_error_log/(2*pi))
xlabel('Time [s]')
ylabel('Error [Hz]')

figure('Name','PLL Error')
plot(time, pll_error_log/(2*pi))
xlabel('Time [s]')
ylabel('Error [cycles]')

figure('Name','Replica Frequency')
hold on
plot(time, fll_replica_freq_log, 'DisplayName','FLL')
plot(time, pll_replica_freq_log, 'DisplayName','PLL')
legend()
xlabel('Time [s]')
ylabel('NCO Frequency [Hz]')

function error = freq_disc(ip, ip_last, qp, qp_last, T)
    cross = ip_last * qp - ip * qp_last;
    dot = ip_last * ip + qp_last * qp;
    error = atan2(cross, dot) / (T); % scaling this is weird
end

function error = phase_disc(ip, qp)
    error = atan(qp/ip);
end

function [freqOffset, fllWNew, fllXNew] = freq_lock_loop(error, noise_bw, fllW, fllX, T)

    a2 = sqrt(2);
    omega = (4 * a2 * noise_bw) / (1 + a2^2);

    fllWNew = fllW + T * omega^2 * error;
    fllXNew = fllX + T * (0.5 * fllWNew + a2 * error);
    freqOffset = 0.5 * fllXNew;

end

function [freqOffset, pllWNew, pllXNew] = phase_lock_loop(pllError, noise_bw, pllW, pllX, T)
    a2 = sqrt(2);
    omega = (4 * a2 * noise_bw) / (1 + a2^2);

    pllWNew = pllW + pllError * omega^2 * T;
    freqOffset = 0.5 * (pllW + pllWNew) + a2 * omega * pllError;
    pllXNew = pllX;
end

function [replica, rem_phase] = create_carr_replica(carr_freq, samp_freq, rem_phase, code_period_size)

    trigTerm = (2 * pi * carr_freq * (1 / samp_freq) * (0:code_period_size)) + rem_phase;
    rem_phase = rem(trigTerm(end), 2*pi);

    replica = exp(-1j*trigTerm(1:code_period_size));
end