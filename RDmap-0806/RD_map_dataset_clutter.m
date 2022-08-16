%
% RD_map_dataset_clutter
%   patch. June 22, 2022
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
H = 2;                   % number of targets
% SNR = [0:5:10];        % dB (need to -30)
SNR = 6;                 % dB (plus 24dB 16*16, 30dB 32*32, 32dB 40*40)
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
% 
% Data matrix construction
% 
for g = 1: length(SNR)
    SNR_g = SNR(g);
    % 
    % target range setting
    % 
    RDmap_signal = cell(TotalSimulationTime,1); % 產生空細胞陣列
    RDmap_noise = cell(TotalSimulationTime,1);
    parameter_setting_record = [];
    for TimeIndex = 1:TotalSimulationTime
        tempMap = zeros(N,M);
        tempMap(1,1) = 1;
        while sum(tempMap(1,:))>0 || sum(tempMap(N,:))>0 || sum(tempMap(:,1))>0 || sum(tempMap(:,M))>0 || sum(sum(tempMap))<H
            fprintf('%d\n', sum(tempMap(1,:)~=0)>0 || sum(tempMap(N,:)~=0)>0 || sum(tempMap(:,1)~=0)>0 || sum(tempMap(:,M)~=0)>0);   
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
            F_Tx = qammod(round(rand(N,M)*4+0.5,0)-1,4)/sqrt(2);
            % 
            % channel effect
            % 
            F_Channel = cell(H,1);
            tempMap = zeros(N,M);
            target_Map(:,:,TimeIndex) = zeros(N,M); % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            for index_target = 1:H
                NuiPhase = rand*2*pi;    % unpredictable phase difference between sources
                range=exp(-1j*2*pi*RoundTripTime(Range(index_target))*SubcarrierSpacing*[1:N].');
                doppler=exp(1j*2*pi*PeriodOFDMsymbol_whole*FreqDopp(Vdop(index_target))*[1:M]);
                F_Channel{index_target} =  F_Tx.*(range*doppler)*exp(1j*NuiPhase);
                RD_map_single_pure_target = abs(fft2(F_Channel{index_target}./F_Tx,N,M)); % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                [nn,mm]=find(RD_map_single_pure_target == max(max(RD_map_single_pure_target)));
                tempMap(nn,mm) = 1; % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
        F_Rx_phase = F_Rx_n./F_Tx;
        Z_processed = Z./F_Tx;
        F_Rx_phase_signal_only = F_Rx./F_Tx;
        RDmap_signal{TimeIndex} = fft2(F_Rx_phase_signal_only,N,M)/sqrt(N)/sqrt(M);   
        RDmap_noise{TimeIndex} = fft2(Z_processed,N,M)/sqrt(N)/sqrt(M);
        % 
        % 沒有 clutter_size_range, clutter_size_Doppler, clutter_power 的資料
        %
        clutter_size_range   = 1;
        clutter_size_Doppler = 1;
        clutter_power = 2;
        
        clutter_range_start = randi(N-clutter_size_range+1); % 沒有variable 'clutter_size_range'.
        clutter_range_end = clutter_range_start+clutter_size_range-1;
        
        clutter_Doppler_start = randi(M-clutter_size_Doppler+1); % 沒有variable 'clutter_size_Doppler'.
        clutter_Doppler_end = clutter_Doppler_start+clutter_size_Doppler-1;
        
        clutter = zeros(N,M);
        clutter(clutter_range_start:clutter_range_end,clutter_Doppler_start:clutter_Doppler_end) = sqrt(clutter_power/2)+1j*sqrt(clutter_power/2);
        % 
        RDmap_full{TimeIndex} = RDmap_signal{TimeIndex} + RDmap_noise{TimeIndex} + clutter;
        RDmap_full_noclutter{TimeIndex} = RDmap_signal{TimeIndex} + RDmap_noise{TimeIndex};
        RDmap_DLlabel_noiseclutter{TimeIndex} = RDmap_noise{TimeIndex} + clutter;
    end
    %
    % vectorization
    %
    for TimeIndex = 1:TotalSimulationTime
        RDmap_input(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_full{TimeIndex}),1,N*M);
        RDmap_input_noclutter(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_full_noclutter{TimeIndex}),1,N*M);
        RDmap_label(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_noise{TimeIndex}),1,N*M); % noise
        RDmap_DLlabel_clutter(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_DLlabel_noiseclutter{TimeIndex}),1,N*M);
        RDmap_label_true_target(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(target_Map(:,:,TimeIndex),1,N*M); % target
    end
end
for i = 1 :  length(SNR)*TotalSimulationTime
    RDmap_input_raw(i,:) = RDmap_input(i,:).^2;
    RDmap_input_raw_noclutter(i,:) = RDmap_input_noclutter(i,:).^2;
    RDmap_label_raw(i,:) = RDmap_label(i,:).^2; % noise
    RDmap_DLlabel_clutter_raw(i,:) = RDmap_DLlabel_clutter(i,:).^2; % noise+clutter
end
% 
% reshape N*M
% 
for ii = 1:size(RDmap_label_true_target,1)
    RD_map_label(:,:,ii) = reshape(RDmap_label_true_target(ii,:),[N,M]); % 只有target
end
for ii = 1:size(RDmap_input_raw,1)% size(A,1):A有幾列 = 模擬次數
    RD_map_clutter(:,:,ii) = reshape(RDmap_input_raw(ii,:),[N,M]); % target+noise+clutter
end
for ii = 1:size(RDmap_input_raw,1)% size(A,1):A有幾列 = 模擬次數
    RD_map_noclutter(:,:,ii) = reshape(RDmap_input_raw_noclutter(ii,:),[N,M]); % target+noise
end
% 
% % % % % 
nmax=zeros(H,size(RDmap_input_raw_noclutter,1));
mmax=zeros(H,size(RDmap_input_raw_noclutter,1));
RD_map_noclutter_2=RD_map_noclutter;
RD_map_noclutter_max=zeros(N,M,TotalSimulationTime);
for jj = 1:size(RDmap_input_raw_noclutter,1)
    for h = 1:H
        [nmax(h,jj),mmax(h,jj)] = find(RD_map_noclutter_2(:,:,jj)==max(max(RD_map_noclutter_2(:,:,jj))));
        RD_map_noclutter_2(nmax(h,jj),mmax(h,jj),jj) = 0;
        RD_map_noclutter_max(nmax(h,jj),mmax(h,jj),jj)=1;
    end   
    [nmaxx(:,jj),mmaxx(:,jj)] = find(RD_map_noclutter_max(:,:,jj)==1);
    [tn(:,jj),tm(:,jj)] = find(RD_map_label(:,:,jj)==1);
    if sum(abs(nmaxx(:,jj)-tn(:,jj))==1)>0 || sum(abs(mmaxx(:,jj)-tm(:,jj))==1)>0
        
        [rdmap_clutter,rdmap_noclutter,label,noise]=newrdmap3(H,SNR);
        RDmap_label_raw(jj,:)=noise;
        RD_map_noclutter(:,:,jj)=rdmap_noclutter;
        temp_rd=RD_map_clutter(:,:,jj);
        RD_map_clutter(:,:,jj)=rdmap_clutter;
        temp_la=RD_map_label(:,:,jj);
        RD_map_label(:,:,jj)=label;
        fprintf('%d,',jj)
    end     
end
% % % % % 
% 
% truncated (optional)
% 
truncated_threshold = 10;
for iii = 1 : TotalSimulationTime
    RDmap_truncated(:,:,iii)=RD_map_clutter(:,:,iii);
    RDmap_input_raw_truncated(iii,:) = reshape(RDmap_truncated(:,:,iii),[],N*M);
    RDmap_input_raw_truncated(iii,find(RDmap_input_raw_truncated(iii,:)>=truncated_threshold)) = truncated_threshold;
%     [xx,yy]=find(RD_map_clutter(:,:,iii)>=truncated_threshold);
%     RD_map_truncated(xx,yy,iii) = truncated_threshold;
    RDmap_DLlabel_clutter_raw(iii,find(RDmap_DLlabel_clutter_raw(iii,:)>=truncated_threshold)) = truncated_threshold;
end  
for ii = 1:size(RDmap_input_raw_truncated,1)
    RD_map_truncated(:,:,ii) = reshape(RDmap_input_raw_truncated(ii,:),[N,M]); 
end
% 
% Dynamic range compression
% 
% % % % % 
for ii = 1:size(RD_map_clutter,3) % 模擬次數
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
% % % % % 
% 
% RD_map
%

% figure(1) % 只有target
% see1 = RD_map_label(:,:,1);
% mesh(see1,'edgecolor','r');

% figure(2) % 只有noise
% see2 = reshape(RDmap_label_raw(1,:),N,M);
% mesh(see2,'edgecolor','r');
% zlim([0,10])

% figure(3) % 沒有clutter
% see3 = RD_map_noclutter(:,:,1);
% mesh(see3,'edgecolor','r');

% figure(4) % DL label noise+clutter
% see7 = reshape(RDmap_DLlabel_clutter_raw(1,:),N,M);
% mesh(see7,'edgecolor','r');

% figure(5) % 有clutter
% see6 = RD_map_clutter(:,:,1);
% mesh(see6,'edgecolor','r');

% figure(6) % 經過truncated
% see4 = reshape(RDmap_input_raw_truncated(1,:),N,M);
% mesh(see4,'edgecolor','r');
% title('truncated')

figure(7) % Dynamic range compression
see5 = reshape(RDmap_input_raw_softknee(1,:),N,M);
mesh(see5,'edgecolor','r');
title('Dynamic range compression: softknee');
figure(8);
imagesc(see5);

% 
% validation
% 
% save(['validation_noT_noclutter_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_noclutter'); % for CFAR
% save(['D:\YU\my_yolo3\DLCFAR\validation_noT_truncated_CNR' num2str(CNR) '_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_truncated'); % for DL-CFAR input
% save(['D:\YU\my_yolo3\DLCFAR\validation_noT_label_noise_CNR' num2str(CNR) '_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_DLlabel_clutter_raw'); % for DL-CFAR label(noise)

% for jj=1:size(RD_map_truncated,3)
%     fileName = ['D:\YU\my_yolo3\VOC2007\JPEGImages_val_trun\validation_noT_truncated_clutter' num2str(CNR) '_H' num2str(H) '_SNR' num2str(SNR) '_f' num2str(jj) '.mat'];
%     RD_map_yoloinput=RD_map_truncated(:,:,jj);
%     save(fileName,'RD_map_yoloinput'); % for YOLO-CFAR input
% end
% save(['validation_noT_label_target_trun_clutter' num2str(CNR) '_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_label'); % for YOLO-CFAR label(target)

% 
% training
% 
% save(['training_noT_noclutter_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_noclutter'); % for CFAR
% save(['D:\YU\my_yolo3\DLCFAR\training_noT_truncated_CNR' num2str(CNR) '_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_truncated'); % for DL-CFAR input
% save(['D:\YU\my_yolo3\DLCFAR\training_noT_label_noise_CNR' num2str(CNR) '_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_DLlabel_clutter_raw'); % for DL-CFAR label(noise)
% % 
% for jj=1:size(RD_map_truncated,3)
%     fileName = ['D:\YU\my_yolo3\VOC2007\JPEGImages_trun\training_noT_truncated_clutter' num2str(CNR) '_H' num2str(H) '_SNR' num2str(SNR) '_f' num2str(jj) '.mat'];
%     RD_map_yoloinput=RD_map_truncated(:,:,jj);
%     save(fileName,'RD_map_yoloinput'); % for YOLO-CFAR input
% end
% save(['training_noT_label_target_trun_clutter' num2str(CNR) '_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_label'); % for YOLO-CFAR label(target)

% 
% testing (H1_30dB被覆蓋)
% 
% save(['testing_clutter_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_clutter'); % for CFAR
% save(['D:\YU\my_yolo3\DLCFAR_clutter\testing_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_truncated'); % for DL-CFAR input
% save(['testing_label_target_clutter_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_label'); % label(target)


% save(['testing_noT_clutter_CNR20_20_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_clutter'); % for CFAR
% save(['D:\YU\my_yolo3\DLCFAR_clutter\testing_noT_truncated_CNR20_20_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_truncated'); % for DL-CFAR input
% save(['D:\YU\my_yolo3\DLCFAR_clutter\testing_noT_label_noise_CNR20_20_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_label_raw'); % for DL-CFAR label(noise)
% for jj=1:size(RD_map_softknee,3)
%     fileName = ['D:\YU\my_yolo3\mAP-master\input\images-optional(20dB_H2_200000_clutter_CNR20)\f' num2str(jj) '_H' num2str(H) '_SNR' num2str(SNR) '.mat'];
%     RD_map_yoloinput=RD_map_softknee(:,:,jj);
%     save(fileName,'RD_map_yoloinput'); % for YOLO-CFAR input
% end
% 

% % YOLO-CFAR label(target)
% save(['testing_noT_label_target_clutter_CNR20_20_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_label');  

% % DNN input
% save(['testing_trainDNN_softknee_clutter_CNR20_20_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_softknee'); 
% save para H SNR N M;
% 

