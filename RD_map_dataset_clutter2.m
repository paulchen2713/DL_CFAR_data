%
% RD_map_dataset_clutter
% 
clear;
clc;
%
% parameter_setting_record = [range1; range2; Vdop1; Vdop2];
% 
% RDmap_signal
% RDmap_noise
%
% Simulation settings
% 
TotalSimulationTime = 1; % number of RD maps
CarrierFreq = 78*10^9;   % carrier frequency 78 GHz
H = 1;                   % number of targets
% SNR = [0:5:10];        % dB (need to -30)
SNR = 0;                 % dB (plus 24dB 16*16, 30dB 32*32, 32dB 40*40)
%
% MIMO parameter settings
%
numTx = 1;
numRx = 1;
% 
% Frame-related parameters
% 
% RD map with size = N*M
N = 4;       % number of subcarrier
M = N;        % number of OFDM symbol
Nfft  = N;    % number of FFT points in frequency domain
frame = 1;    % number of frame used
Mfft  = M;    % number of FFT points in time domain
%
% OFDM-related parameters
%
BW = 1*10^9;  % 1 GHz
SubcarrierSpacing = BW/N;
PeriodOFDMsymbol = 1/SubcarrierSpacing;
CP = 49*PeriodOFDMsymbol;
PeriodOFDMsymbol_whole = PeriodOFDMsymbol + CP;
PeriodFrame = PeriodOFDMsymbol_whole * M;
% 
% Common parameter definition
% 
c = 299792458;    % speed of light
FreqDopp = @(RelVelocity) 2*RelVelocity*CarrierFreq/c; % Doppler shift function
RoundTripTime = @(d) 2*d/c;                            % Round Trip Time (RTT) function
%
% Specifications
%
d_unamb = c/2/SubcarrierSpacing;                       % unambiguity range
v_unamb = c/2/CarrierFreq/PeriodOFDMsymbol_whole;      % unambiguity velocity
d_rel = d_unamb/N; % search resolution of range
v_rel = v_unamb/M; % search resolution of velocity
% 
% Data matrix construction
% 

% 
% target range setting
% 
RDmap_signal = cell(TotalSimulationTime,1); % 玻ネ灿璏皚
RDmap_noise = cell(TotalSimulationTime,1);
parameter_setting_record = [];

for h = 1:H
    Range(h,1) = rand*d_unamb;
end
% target Doppler velocity setting
for h = 1:H
    Vdop(h,1) = (2*rand-1)*v_unamb;
end
% 
% target DoA setting
% 
for h = 1:H
    DoA(h,1) = (2*rand-1)*60; % FOV = 120 degrees
end
parameter_setting_record = [parameter_setting_record [Range;Vdop]]; 

% 
% transmitted signal (QPSK symbols)
% 
F_Tx = qammod(round(rand(N,M)*4 + 0.5, 0) - 1, 4) / sqrt(2);
% 
% channel effect
% 
F_Channel = cell(H,1);
for index_target = 1:H
    NuiPhase = rand*2*pi;    % unpredictable phase difference between sources
    range=exp(-1j*2*pi*RoundTripTime(Range(index_target))*SubcarrierSpacing*(1:N).');
    doppler=exp(1j*2*pi*PeriodOFDMsymbol_whole*FreqDopp(Vdop(index_target))*(1:N));
    F_Channel{index_target} =  F_Tx .* (range*doppler)*exp(1j*NuiPhase);

    fft2_input = F_Channel{index_target}./F_Tx;
    fft2_result = fft2(fft2_input, N, M);

    RD_map_single_pure_target = abs(fft2_result); % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    [nn,mm] = find(RD_map_single_pure_target == max(max(RD_map_single_pure_target)));
end 
% 
% received signal
% 
F_Rx = cell(1,1);
F_Rx_phase = cell(1,1);    
F_Rx = zeros(N,M);        
P_noise = 0.5;
P_signal = P_noise*(10^(SNR/10));
for index_target = 1:H
    F_Rx = F_Rx + sqrt(P_signal/2)*F_Channel{index_target};
end
% 
% add sptially white noise
% 
P_Rx = mean(mean(F_Rx.*conj(F_Rx)));
Z = sqrt(P_noise/2)*(randn(N,M)+1j*randn(N,M));
F_Rx_n = F_Rx + Z;
F_Rx_phase = F_Rx_n./F_Tx;
Z_processed = Z./F_Tx;
F_Rx_phase_signal_only = F_Rx ./ F_Tx;
RDmap_signal = fft2(F_Rx_phase_signal_only, N, M)/sqrt(N)/sqrt(M);   
RDmap_noise = fft2(Z_processed,N,M)/sqrt(N)/sqrt(M);
% 
% ⊿Τ clutter_size_range, clutter_size_Doppler, clutter_power 戈
%
clutter_size_range   = 1;
clutter_size_Doppler = 1;
clutter_power = 2;

clutter_range_start = randi(N-clutter_size_range+1); % ⊿Τvariable 'clutter_size_range'.
clutter_range_end = clutter_range_start+clutter_size_range-1;

clutter_Doppler_start = randi(M-clutter_size_Doppler+1); % ⊿Τvariable 'clutter_size_Doppler'.
clutter_Doppler_end = clutter_Doppler_start+clutter_size_Doppler-1;

clutter = zeros(N,M);
clutter(clutter_range_start:clutter_range_end,clutter_Doppler_start:clutter_Doppler_end) = sqrt(clutter_power/2)+1j*sqrt(clutter_power/2);
% 
RDmap_full = RDmap_signal + RDmap_noise + clutter;
RDmap_full_noclutter = RDmap_signal + RDmap_noise;

%
% vectorization
%
RDmap_input = reshape(abs(RDmap_full),1,N*M);
RDmap_input_noclutter = reshape(abs(RDmap_full_noclutter),1,N*M);

RDmap_input_raw(:) = RDmap_input(:).^2;
RDmap_input_raw_noclutter(:) = RDmap_input_noclutter(:).^2;
% 
% reshape N*M
% 
for ii = 1:size(RDmap_input_raw,1)% size(A,1):AΤ碭 = 家览Ω计
    RD_map_clutter(:,:,ii) = reshape(RDmap_input_raw(ii,:),[N,M]); % target+noise+clutter
end
for ii = 1:size(RDmap_input_raw,1)% size(A,1):AΤ碭 = 家览Ω计
    RD_map_noclutter(:,:,ii) = reshape(RDmap_input_raw_noclutter(ii,:),[N,M]); % target+noise
end


% % % % % 
% 
% truncated (optional)
% 
truncated_threshold = 10;

figure(8);
imagesc(RD_map_clutter);
figure(9);
mesh(RD_map_clutter,'edgecolor','r');

RDmap_input_raw_truncated(1,:) = reshape(RD_map_clutter(:,:,1),[],N*M);
RDmap_input_raw_truncated(1, find(RDmap_input_raw_truncated(1,:) >= truncated_threshold) ) = truncated_threshold;

for ii = 1:size(RDmap_input_raw_truncated,1)
    RD_map_truncated(:,:,ii) = reshape(RDmap_input_raw_truncated(ii,:),[N,M]); 
end

figure(10);
imagesc(RD_map_truncated);
figure(11);
mesh(RD_map_truncated,'edgecolor','r');

% 
% Dynamic range compression
% 
% % % % % 
for ii = 1:size(RD_map_clutter,3)% 家览Ω计
    RDmap_input_raw_clutter2(ii,:) = reshape(RD_map_clutter(:,:,ii),[],N*M);
end
threshhold_noise = 4;
knee_stop = 20;
Maximum = 30;
R = 88;
T = 12;
W = 16;
RDmap_input_raw_softknee = RDmap_input_raw_clutter2;
RDmap_input_raw_log = RDmap_input_raw_clutter2;
for i = 1 : TotalSimulationTime
    for j = 1: N*M
        if threshhold_noise < RDmap_input_raw_clutter2(i,j) && RDmap_input_raw_clutter2(i,j)  < knee_stop
            RDmap_input_raw_softknee(i,j) = RDmap_input_raw_clutter2(i,j)+(1/R-1)*((RDmap_input_raw_clutter2(i,j) -T+W/2).^2)/(2*W);
        elseif RDmap_input_raw_clutter2(i,j)  > knee_stop
            RDmap_input_raw_softknee(i,j) = T + (RDmap_input_raw_clutter2(i,j) - T)/R;
        end
    end
end
% 
for ii = 1:size(RDmap_input_raw_softknee,1)
    RD_map_softknee(:,:,ii) = reshape(RDmap_input_raw_softknee(ii,:),[N,M]); % Dynamic range compression: softknee
end

figure(12);
imagesc(RD_map_softknee);
figure(13);
mesh(RD_map_softknee,'edgecolor','r');
