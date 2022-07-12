%
% save raw RDmaps
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
SNR = 0:5:10;            % dB (need to -30)
% SNR = 6;                 % dB (plus 24dB 16*16, 30dB 32*32, 32dB 40*40)
%
% MIMO parameter settings
%
numTx = 1;
numRx = 1;
% 
% Frame-related parameters
% 
% RD map with size = N*M
N = 16;       % number of subcarrier
M = 16;       % number of OFDM symbol
Nfft  = N;    % number of FFT points in frequency domain
frame = 1;    % number of frame used
Mfft  = M;    % number of FFT points in time domain
%
% OFDM-related parameters
%
BW = 1*10^9;
SubcarrierSpacing = BW/N;
PeriodOFDMsymbol = 1/SubcarrierSpacing;
CP = 49*PeriodOFDMsymbol;
PeriodOFDMsymbol_whole = PeriodOFDMsymbol + CP;
PeriodFrame = PeriodOFDMsymbol_whole * M;
% 
% Common parameter definition
% 
c = 3*10^8;    % speed of light
FreqDopp = @(RelVelocity) 2*RelVelocity*CarrierFreq/c; % Doppler shift function
RoundTripTime = @(d) 2*d/c;                            % Round Trip Time (RTT) function
%
% Specifications
%
d_unamb = c/2/SubcarrierSpacing;                       % unambiguity range
v_unamb = c/2/CarrierFreq/PeriodOFDMsymbol_whole;      % unambiguity velocity
d_rel = d_unamb/N; % search resolution of range
v_rel = v_unamb/M; % search resolution of velocity

max_iter = 1;
for iter = 1:max_iter
    fprintf("iter: %d\n", iter);
    % 
    % Data matrix construction
    %
    for g = 1: length(SNR)
        SNR_g = SNR(g);
        % 
        % target range setting
        % 
        RDmap_signal = cell(TotalSimulationTime,1); % 玻ネ灿璏皚
        RDmap_noise  = cell(TotalSimulationTime,1);
        parameter_setting_record = [];
        for TimeIndex = 1:TotalSimulationTime
            %
            tempMap = zeros(N,M);
            tempMap(1,1) = 1;
            while sum(tempMap(1,:))>0 || sum(tempMap(N,:))>0 || sum(tempMap(:,1))>0 || sum(tempMap(:,M))>0 || sum(sum(tempMap))<H
                fprintf('%d\n', sum(tempMap(1,:)~=0)>0 || sum(tempMap(N,:)~=0)>0 || sum(tempMap(:,1)~=0)>0 || sum(tempMap(:,M)~=0)>0);   
                for h = 1:H
                    Range(h,1) = rand*d_unamb;
                end
                % 
                % target Doppler velocity setting
                % 
                for h = 1:H
                    Vdop(h,1) = (2*rand-1)*v_unamb;
                end
                % 
                % target DoA setting
                % 
                for h = 1:H
                    DoA(h,1) = (2*rand-1)*60; % FOV = 120 degrees
                end
                parameter_setting_record = [parameter_setting_record [Range; Vdop]]; 
                % 
                % transmitted signal (QPSK symbols)
                % 
                F_Tx = qammod(round(rand(N,M)*4+0.5,0)-1,4)/sqrt(2);  
                % 
                % channel effect
                % 
                F_Channel = cell(H,1);
                tempMap = zeros(N,M);
                target_Map(:,:,TimeIndex) = zeros(N,M); % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                for index_target = 1:H
                    % unpredictable phase difference between sources
                    NuiPhase = rand*2*pi;
                    range    = exp(-1j*2*pi*RoundTripTime(Range(index_target))*SubcarrierSpacing*[1:N].');
                    doppler  = exp(1j*2*pi*PeriodOFDMsymbol_whole*FreqDopp(Vdop(index_target))*[1:M]);
                    %
                    F_Channel{index_target} =  F_Tx.*(range*doppler)*exp(1j*NuiPhase);
                    RD_map_single_pure_target = abs(fft2(F_Channel{index_target}./F_Tx,N,M)); % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                    % 
                    [nn, mm] = find(RD_map_single_pure_target == max(max(RD_map_single_pure_target)));
                    tempMap(nn, mm) = 1; % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                end       
            end
            target_Map(:,:,TimeIndex) = tempMap; % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % 
            % received signal
            % 
            F_Rx = cell(1,1);
            F_Rx_phase = cell(1,1);
            F_Rx = zeros(N,M);        
            P_noise = 0.5;
            P_signal = P_noise*(10^(SNR_g/10));
            for index_target = 1:H
                F_Rx = F_Rx + sqrt(P_signal/2)*F_Channel{index_target};
            end
            % 
            % add sptially white noise
            % 
            P_Rx = mean(mean(F_Rx.*conj(F_Rx)));
            Z = sqrt(P_noise/2)*(randn(N,M)+1j*randn(N,M));
            F_Rx_n = F_Rx + Z;
            F_Rx_phase  = F_Rx_n./F_Tx;
            Z_processed = Z./F_Tx;
            F_Rx_phase_signal_only  = F_Rx./F_Tx;
            RDmap_signal{TimeIndex} = fft2(F_Rx_phase_signal_only,N,M)/sqrt(N)/sqrt(M);
            RDmap_noise{TimeIndex}  = fft2(Z_processed,N,M)/sqrt(N)/sqrt(M);
            RDmap_full_noclutter{TimeIndex} = RDmap_signal{TimeIndex} + RDmap_noise{TimeIndex};
        end
        % 
        % vectorization
        % 
        for TimeIndex = 1:TotalSimulationTime
            RDmap_input_noclutter(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_full_noclutter{TimeIndex}),1,N*M);
            RDmap_label(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_noise{TimeIndex}),1,N*M); % noise
            RDmap_label_true_target(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(target_Map(:,:,TimeIndex),1,N*M); % target
        end
    end
    %
    % truncated (optional)
    %
    truncated_threshold = 10;
    for i = 1 :  length(SNR)*TotalSimulationTime
        RDmap_input_raw_noclutter(i,:) = RDmap_input_noclutter(i,:).^2;
        RDmap_label_raw(i,:) = RDmap_label(i,:).^2;
        RDmap_input_raw_truncated(i,:) = RDmap_input_raw_noclutter(i, :);
        RDmap_input_raw_truncated(i, find(RDmap_input_raw_noclutter(i, :) >= truncated_threshold)) = truncated_threshold;
    end
    %
    % reshape N*M
    %
    for ii = 1:size(RDmap_label_true_target,1)
        RD_map_label(:,:,ii) = reshape(RDmap_label_true_target(ii,:), [N,M]); % Τtarget
    end
    for ii = 1:size(RDmap_input_raw_noclutter,1) % size(A,1):AΤ碭 = 家览Ω计
        RD_map_noclutter(:,:,ii) = reshape(RDmap_input_raw_noclutter(ii, :), [N,M]); % target+noise
    end
    %


    %
    % Dynamic Range Compression (DRC)
    %
    for ii = 1:size(RD_map_noclutter,3) % 家览Ω计
        RDmap_input_raw_noclutter2(ii,:) = reshape(RD_map_noclutter(:,:,ii),[],N*M);
    end
    threshhold_noise = 4;
    knee_stop = 20;
    Maximum = 30;
    R = 88;
    T = 12;
    W = 16;
    RDmap_input_raw_softknee = RDmap_input_raw_noclutter2;
    RDmap_input_raw_log = RDmap_input_raw_noclutter2;
    for i = 1 : TotalSimulationTime
        for j = 1: N*M
            if threshhold_noise < RDmap_input_raw_noclutter2(i,j) && RDmap_input_raw_noclutter2(i,j)  < knee_stop
                RDmap_input_raw_softknee(i,j) = RDmap_input_raw_noclutter2(i,j)+(1/R-1)*((RDmap_input_raw_noclutter2(i,j) -T+W/2).^2)/(2*W);
            elseif RDmap_input_raw_noclutter2(i,j)  > knee_stop
                RDmap_input_raw_softknee(i,j) = T + (RDmap_input_raw_noclutter2(i,j) - T)/R;
            end
        end
    end
    for i = 1 : TotalSimulationTime
        for j = 1: N*M
            RDmap_input_raw_log(i,j) = log(RDmap_input_raw_noclutter2(i,j)+1)*30/7;
        end
    end
    % 
    % RD_map_softknee = zeros(16, 16, 3);
    % RD_map_log = zeros(16, 16, 3);
    for ii = 1:size(RDmap_input_raw_softknee,1)
        RD_map_softknee(:,:,ii) = reshape(RDmap_input_raw_softknee(ii,:),[N,M]); % Dynamic range compression: softknee
    end
    for ii = 1:size(RDmap_input_raw_log,1)
        RD_map_log(:,:,ii) = reshape(RDmap_input_raw_log(ii,:),[N,M]); % Dynamic range compression: log
    end
    %
    % RD_map
    %
    % figure(1) % Τtarget
    % see1 = RD_map_label(:,:,1);
    % mesh(see1,'edgecolor','r');
    % 
    % figure(2) % Τnoise
    % see2 = reshape(RDmap_label_raw(1,:),N,M);
    % mesh(see2,'edgecolor','r');
    % zlim([0,10])
    % 
    % figure(3) % ⊿Τclutter
    % see3 = reshape(RDmap_input_raw_noclutter2(1,:),N,M);
    % mesh(see3,'edgecolor','r');
    % 
    % figure(10)
    % imagesc(see3);
    % 
    % figure(4) % 竒筁truncated
    % see4 = RD_map_truncated(:,:,1);
    % mesh(see4,'edgecolor','r');
    % title('Truncated')
    % 
    % figure(5) % Dynamic range compression
    see5 = reshape(RDmap_input_raw_softknee(1, :), N, M);
    % mesh(see5,'edgecolor','r');
    % title('Dynamic range compression');
%     figure(iter)
    imagesc(see5);
    
    set(gca,'XTick',[]) % Remove the ticks in the x axis
    set(gca,'YTick',[]) % Remove the ticks in the y axis
    set(gca,'Position', [0 0 1 1]) % Make the axes occupy the hole figure
    
    % % saveas(figureHandle,'filename','format') 
    % gcf means "get current figure"
    cur_filename = [num2str(iter),'.jpeg'];
    saveas(gcf, cur_filename, 'jpeg');
end


