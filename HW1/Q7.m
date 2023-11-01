
%% Packet Construction
clear
clc
B = 20e6; % Sampling rate 
N = 64; % FFT Size
L = 16; % Cyclic Prefix
%%step a)
%% Get random vector contains 4176 bits
my_packet =  randi([0,1], 1, 4176);

%% Convert them into BPSK symbols
my_BPSK_symbols = zeros(1,4176);
for i = 1: length(my_packet)
    element = my_packet(i);
    if  element == 0
        my_BPSK_symbols(i) = 1 + 0i;
    else
        my_BPSK_symbols(i) = -1 + 0i;
    end
end
N_data = 48; % number of data
my_BPSK_symbols = complex(my_BPSK_symbols);
N_block = ceil(length(my_BPSK_symbols) / N_data);
my_BPSK_symbols = reshape(my_BPSK_symbols,N_data, N_block);
my_BPSK_symbols = my_BPSK_symbols';
%% Group BPSK symbols into 802.11 OFDM symbols

N_pilots = 4; % number of pilots
pilots = 1+0i; % pilots
pilots_index = [7, 21, 43, 57]; %piltos index in the ofdm dymbol
pilots_index = pilots_index + 1; %piltos index in the ofdm dymbol

N_nulls = 12; % number of nulls
nulls_index = [0, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]; % nulls' index
nulls_index = nulls_index + 1; % nulls' index

data_index = setdiff(1:N, [pilots_index, nulls_index]);

my_OFDM_symbol = complex(zeros(N_block, N));% 87 by 64

% Insert data
% 87 by 64
my_OFDM_symbol(:,data_index) = my_BPSK_symbols; 
% Insert pilot
my_OFDM_symbol(:,pilots_index) = 1;
% Insert null
my_OFDM_symbol(:,nulls_index) = 0;
%% step b) OFDM Modulation

N_cyclic_prefix = 16;% 16 samples for cyclic prefix
% 64-point IFFT
my_modulation_time = ifft(my_OFDM_symbol, N, 2);

% Add a 16 sample cyclic prefix
% 87 by 80
my_modulation_time_prefix = cat(2,my_modulation_time(:,end-N_cyclic_prefix+1:end),my_modulation_time);

%% step c) Add STF and LTF preamables
stf_freq = complex(zeros(1, N)); % From 802.11
stf_freq(1,end-25:end) = sqrt(13/6) * [0, 0, 1+1i, 0, 0, 0, -1-1i, 0, 0, 0, 1+1i, 0, 0, 0, -1-1i, 0, 0, 0, -1-1i, 0, 0, 0, 1+1i, 0, 0, 0];
stf_freq(1, 1:27) = sqrt(13/6) * [0, 0, 0, 0, -1-1i, 0, 0, 0, -1-1j, 0, 0, 0, 1+1i, 0, 0, 0, 1+1i, 0, 0, 0, 1+1i, 0, 0, 0, 1+1i, 0,0];
stf_time = ifft(stf_freq, N,2);
stf_symbol = stf_time(1:16);
stf_all = cat(2, stf_symbol,stf_symbol, stf_symbol,stf_symbol, stf_symbol,stf_symbol, stf_symbol,stf_symbol, stf_symbol,stf_symbol);

ltf_freq = complex(zeros(1,N));
ltf_freq(1,end-25:end) = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1];
ltf_freq(1,1:27) = [0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];
ltf_time = ifft(ltf_freq, N,2);
ltf_symbol = ltf_time;
ltf_all = cat(2, ltf_symbol(end-31: end), ltf_symbol, ltf_symbol);


% Concatencate STF and LTF
prefix_data = reshape(my_modulation_time_prefix', 1, 6960);
%% step d) Plot STF and Power Spectrum Density of the OFDM data symbols
figure('Name','Magnitude of STF');
plot(linspace(1,160,160), abs(stf_all));
xlabel('index');
ylabel('Magnitude');
title('Magnitude of STF');
saveas(gcf, 'Magnitude of STF.png');

figure('Name','Power Density of the OFDM')
plot(linspace(1, N, N),20*log10(abs(fft(prefix_data,64,2))));
xlabel('frequency');
ylabel('Power(dB)');
title('Power Densityof the OFDM');
saveas(gcf, 'Power_Spectrum_Density.png');
%% Packet transmission and channel distortion
data_stf_ltf = cat(2, stf_all, ltf_all,prefix_data);
data_with_idle = cat(2, zeros(1,100), data_stf_ltf); % Add idle period
% Define variables
magnitude_distortion = 10^(-5);
channel_shift = exp(-1i*3*pi/4);
freq_offset = exp(-1i*2*pi*0.00017 * (1:length(data_with_idle)));
channel_noise = 10^(-14);
random_noise = randn(1,length(data_with_idle));
random_noise = random_noise*channel_noise + 0;
% After Distortion
data_after_distortion = data_with_idle *magnitude_distortion * channel_shift.* freq_offset + random_noise;

figure('Name','Magnitude_of_STF_after_distortion_effects');
plot(linspace(1,160,160), abs(data_after_distortion(:, 101:101+160-1)));
xlabel('index');
ylabel('Magnitude');
title('Magnitude of STF after distortion effects');
saveas(gcf, 'Magnitude_of_STF_after_distortion_effects.png');
%% Packet Detection
len_Detection = 500; % Known stf is in 101-260
correlation = zeros(1, len_Detection);
energy = zeros(1,len_Detection);
for m=1:len_Detection
    % Correlation = current 16 samples and next 16samples
    correlation(1,m) = dot(data_after_distortion(:,m:m+16-1), data_after_distortion(:,m+16 :m+16 + 15));
    % Energy: E = data^2
    energy(1, m) = dot(data_after_distortion(:,m:m+16-1), data_after_distortion(:,m:m+16-1));
end
figure('Name','Packet Detection');
plot(linspace(1,len_Detection,len_Detection), abs(correlation),'b-');
hold on 
plot(linspace(1,len_Detection,len_Detection), abs(energy),'r--');
xlabel('index');
ylabel('Magnitude');
legend('correlation', 'energy');
title('Packet Detection');
saveas(gcf, 'Packet_Detection.png');
% Find where packet is detected: when they are roughly the same
stf_find_detect = find(abs(correlation) > 0.95*abs(energy)) + 31
% 101-260

%% Packet Synchronization
% Cross-correlation, detected if the know sequence has occurred
cross_correlation = zeros(1,500);
% channel_distortion_synchronization = zeros(1,500);
for m = 1: 500
     cross_correlation(1, m) = dot(data_after_distortion(:,m:m+16-1), stf_all(:,1:16));
     % channel_distortion_synchronization(1,m) = dot(magnitude_distortion * channel_shift.* freq_offset(:,m:m+16-1), magnitude_distortion * channel_shift.* freq_offset(:,m:m+16-1));
end
figure('Name','Packet Synchronization, cross-correlation');
plot(linspace(1,500,500), abs(cross_correlation),'b-');
% hold on 
% plot(linspace(1,len_Detection,len_Detection), abs(channel_distortion_synchronization),'r--');
xlabel('index');
ylabel('Magnitude');
legend('cross-correlation');
title('Packet Synchronization, cross-correlation');
saveas(gcf, 'Packet-Synchronization_cross-correlation.png');
% Find where stf started
str_find_synch = find(abs(cross_correlation) >= 0.95*max(abs(cross_correlation)))
% 101   117   133   149   165   181   197   213   229   245
%%  Channel estimation and packet decoding
delta_freq = zeros(1, 160);
for i = 261 : 261+159-64 % 261(100 + 160) where ltf start
    y_n = data_after_distortion(i);
    y_n_N = data_after_distortion(i + N);
    delta_freq(1, i) = abs((imag(y_n_N / y_n)) / (2*pi *64));
end
delta_freq = delta_freq(1,261:end);
delta_freq_avg = sum(delta_freq)/96 %1.6987e-04
fprintf("Frequency_offset: %d\n", delta_freq_avg); % print the frequency offset
data_after_compensation = data_after_distortion .* exp(1i* 2 * pi* delta_freq_avg * linspace(1,length(data_after_distortion), length(data_after_distortion))); %exp(-j2piif)

% Channel gain estimation

% For LTF, received signal on subcarrier k is Y(k) = FFT(LTF1) or FFT(LTF2)
Y = fft(data_after_compensation(261 + 32:261 + 32 + 63), 64, 2);
fprintf("Distortion:\n"); % print the distortion
H = Y ./ ltf_freq

%% Decode
% starting time of OFDM data symbol is 100 + 320
OFDM_data_symbol = data_after_compensation(100 + 320 + 1: end); % 1 by 6960

% Remove cyclic prefix
OFDM_data_symbol_reshape = (reshape(OFDM_data_symbol, 80, 87))'; % 87 by 80

OFDM_data_symbol_noprefix = OFDM_data_symbol_reshape(:,17:end); % 87 by 64

OFDM_data_symbol_noprefix_freq = fft(OFDM_data_symbol_noprefix, 64, 2); % fft 64

data_revert = zeros(87, 64);
H_repeated = repmat(H, 87,1);
data_revert(:, data_index) = OFDM_data_symbol_noprefix_freq(:, data_index) ./ H_repeated(:, data_index); % Y(k) / H(k)

data_revert = data_revert(:, data_index);

% Demap
bpsk_symbol_recover = reshape(data_revert', 4176, 1)';

packet_recovery = zeros(1, 4176);
for i = 1: length(bpsk_symbol_recover)
    if real(bpsk_symbol_recover(i)) < 0
        packet_recovery(:, i) = 1;
    else
        packet_recovery(:, i) = 0;
    end
end

% Verify if they are the same
if packet_recovery == my_packet
    fprintf("I send:\n");
    disp(my_packet);
    fprintf("I received:\n");
    disp(packet_recovery);
    fprintf("Yes you got the same packet\n");
else
    fprintf("Noooooo!!! They are not the same\n"); 
end