%
% RD_map_dataset main program
%   patch. June 22, 2022
%          July 02, 2022
%          July 11, 2022
%
clear;
clc;
%
% parameter_setting_record = [range1; range2; Vdop1; Vdop2];
% 
% RDmap_signal, RDmap_noise
%
% Simulation settings
% 
TotalSimulationTime = 1; % number of RD maps
CarrierFreq = 78*10^9;   % carrier frequency 78 GHz
H = 1;                   % number of targets        但只能設1個target !?
SNR = 0:2:10;            % dB (need to -30)??  0:5:10
% SNR = 6;                 % dB (plus 24dB 16*16, 30dB 32*32, 32dB 40*40)   6dB? 單位很奇怪
%
% MIMO parameter settings 
%
% numTx = 1;
% numRx = 1;
% 
% Frame-related parameters
% 
% RD map with size = N*M
N = 16;       % number of subcarrier
M = 16;       % number of OFDM symbol
% Nfft  = N;    % number of FFT points in frequency domain
% frame = 1;    % number of frame used
% Mfft  = M;    % number of FFT points in time domain
%
% OFDM-related parameters
%
BW = 1*10^9;
SubcarrierSpacing = BW/N;
PeriodOFDMsymbol = 1/SubcarrierSpacing;
CP = 49*PeriodOFDMsymbol; % 49?? CP?? 
PeriodOFDMsymbol_whole = PeriodOFDMsymbol + CP;
PeriodFrame = PeriodOFDMsymbol_whole * M;       % unused
% 
% Common parameter definition
% 
c = physconst('LightSpeed');                          % speed of light == 299792458 ~ 3*10^8 m/s
FreqDopp = @(RelVelocity)2*RelVelocity*CarrierFreq/c; % Doppler shift function
RoundTripTime = @(d)2*d/c;                            % Round Trip Time (RTT) function
%
% Specifications
%
d_unamb = c/2/SubcarrierSpacing;                       % unambiguous range is the maximum range at which a target can be located
v_unamb = c/2/CarrierFreq/PeriodOFDMsymbol_whole;      % unambiguous velocity
% d_resol = d_unamb/N;                                   % search resolution of range
% v_resol = v_unamb/M;                                   % search resolution of velocity
% 
% Data matrix construction
%
for g = 1:length(SNR)
    SNR_g = SNR(g);
    % 
    % target range setting
    % 
    % cell(), called cell array, is a data type with indexed data containers called cells
    % 
    % initialize RDmap_signal and RDmap_noise as empty TotalSimulationTime-by-1 cell arrays
    RDmap_signal = cell(TotalSimulationTime, 1); % cell(1, 1) = 1×1 cell array {0×0 double}
    RDmap_noise  = cell(TotalSimulationTime, 1); 
    % parameter_setting_record = [];
    for TimeIndex = 1:TotalSimulationTime
        %
        tempMap = zeros(N, M); 
        tempMap(1, 1) = 1;     % why initialize first element in tempMap as 1?
        % 
        % 如果 tempMap的第1個或第N個row 的元素和大於0 或是 tempMap的第1個或第N個column 的元素和大於0 
        % 或是 tempMap的總元素和小於 目標數H 任一條件滿足就跳出迴圈, sum(A) returns the sum of the elements
        while sum(tempMap(1, :)) > 0 || sum(tempMap(N, :)) > 0 || sum(tempMap(:, 1)) > 0 || sum(tempMap(:, M)) > 0 || sum(sum(tempMap)) < H
            % fprintf('%d\n', sum(tempMap(1,:)~=0)>0 || sum(tempMap(N,:)~=0)>0 || sum(tempMap(:,1)~=0)>0 || sum(tempMap(:,M)~=0)>0);   
            % 
            Range = zeros(H, 1);
            Vdop  = zeros(H, 1);
            % DoA   = zeros(H, 1);
            %
            for h = 1:H
                % 都卜勒距離 = 隨機值 * 絕對距離, rand ~ U(0, 1)
                Range(h, 1) = rand * d_unamb;
            end
            % 
            % target Doppler velocity setting
            % 
            for h = 1:H
                % 都卜勒速度 = (2*隨機值 - 1) * 絕對速度
                Vdop(h, 1) = (2*rand - 1) * v_unamb;
            end
            % 
            % target DoA setting
            % 
            % for h = 1:H
            %     DoA(h, 1) = (2*rand - 1)*60; % FOV = 120 degrees 
            % end
            % parameter_setting_record = [parameter_setting_record [Range; Vdop]]; 
            % 
            % transmitted signal (QPSK symbols)
            % 
            F_Tx = qammod(round(rand(N, M)*4 + 0.5, 0) - 1, 4) / sqrt(2);  
            % 
            % channel effect
            % 
            F_Channel = cell(H, 1);
            tempMap   = zeros(N, M); % 清空tempMap 
            target_Map(:,:,TimeIndex) = zeros(N,M); % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            for index_target = 1:H
                % unpredictable phase difference between sources
                NuiPhase = rand*2*pi;
                range    = exp(-1j*2*pi*RoundTripTime(Range(index_target))*SubcarrierSpacing*[1:N].');
                doppler  = exp(1j*2*pi*PeriodOFDMsymbol_whole*FreqDopp(Vdop(index_target))*[1:M]);
                %
                F_Channel{index_target} =  F_Tx .* (range*doppler)*exp(1j*NuiPhase);
                RD_map_single_pure_target = abs(fft2(F_Channel{index_target}./F_Tx,N,M)); % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % 
                [nn, mm] = find(RD_map_single_pure_target == max(max(RD_map_single_pure_target)));
                tempMap(nn, mm) = 1; % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            end
        end % end while, 代表 tempMap 有不為零的elemennt, 所以assign給target必不為0
        target_Map(:, :, TimeIndex) = tempMap; % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % 
        % received signal
        % 
        % F_Rx = cell(1,1);
        % F_Rx_phase = cell(1,1);
        F_Rx = zeros(N, M);        
        P_noise = 0.5;      % unint?
        P_signal = P_noise*(10^(SNR_g/10));
        for index_target = 1:H
            F_Rx = F_Rx + sqrt(P_signal/2)*F_Channel{index_target};
        end
        % 
        % add sptially white noise
        % 
        P_Rx = mean(mean(F_Rx.*conj(F_Rx)));
        Z = sqrt(P_noise/2)*(randn(N, M) + 1j*randn(N, M));
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
% truncated (not optional)
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
RD_map_label = zeros(N, M, 3);
for ii = 1:size(RDmap_label_true_target, 1)
    % 
    RD_map_label(:,:,ii) = reshape(RDmap_label_true_target(ii,:), [N, M]); % 只有target
end
for ii = 1:size(RDmap_input_raw_noclutter,1) % size(A,1):A有幾列 = 模擬次數
    RD_map_noclutter(:,:,ii) = reshape(RDmap_input_raw_noclutter(ii, :), [N,M]); % target+noise
end
%
% 
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
    [tn(:,jj), tm(:,jj)] = find(RD_map_label(:,:,jj)==1);
    if sum(abs(nmaxx(:,jj)-tn(:,jj))==1)>0 || sum(abs(mmaxx(:,jj)-tm(:,jj))==1)>0
        % newrdmap2 for single-target?
        % newrdmap3 for multi-target? nope
        [rdmap,label] = newrdmap2(H,SNR); 
        temp_rd=RD_map_noclutter(:,:,jj);
        RD_map_noclutter(:,:,jj)=rdmap;
        temp_la=RD_map_label(:,:,jj);
        RD_map_label(:,:,jj)=label;
        fprintf('%d,',jj)
    end     
end
% 

% % % % % 
% nmax=zeros(H,size(RDmap_input_raw_noclutter,1));
% mmax=zeros(H,size(RDmap_input_raw_noclutter,1));
% RD_map_noclutter_2=RD_map_noclutter;
% RD_map_noclutter_max=zeros(N,M,TotalSimulationTime);
% for jj = 1:size(RDmap_input_raw_noclutter,1)
%     for h = 1:H
%         [nmax(h,jj),mmax(h,jj)] = find(RD_map_noclutter_2(:,:,jj)==max(max(RD_map_noclutter_2(:,:,jj))));
%         RD_map_noclutter_2(nmax(h,jj),mmax(h,jj),jj) = 0;
%         RD_map_noclutter_max(nmax(h,jj),mmax(h,jj),jj)=1;
%     end   
%     [nmaxx(:,jj),mmaxx(:,jj)] = find(RD_map_noclutter_max(:,:,jj)==1);
%     [tn(:,jj),tm(:,jj)] = find(RD_map_label(:,:,jj)==1);
%     if sum(abs(nmaxx(:,jj)-tn(:,jj))==1)>0 || sum(abs(mmaxx(:,jj)-tm(:,jj))==1)>0      
%         fprintf('p%d,',jj)
%     end   
% end 
% % % % % 

% % % % % 
% for ii = 1:size(RDmap_input_raw_noclutter,1)
%     [nmax(ii), mmax(ii)] = find(RD_map_noclutter(:,:,ii) == max(max(RD_map_noclutter(:,:,ii))));
%     [tn(ii), tm(ii)] = find(RD_map_label(:,:,ii)==1);
%     if nmax(ii) ~= tn(ii) || mmax(ii) ~= tm(ii)
%         fprintf('0\n')
%     else
%         fprintf('3\n')
%     end
% end
% % % % % 

%
% truncated (optional)
%
% truncated_threshold = 10;
% for iii = 1 : TotalSimulationTime
%     RD_map_truncated(:,:,iii) = RD_map_noclutter(:,:,iii);
%     [xx,yy]=find(RD_map_noclutter(:,:,iii)>=truncated_threshold);
%     for xxx = 1:length(xx)
%         RD_map_truncated(xx(xxx),yy(xxx),iii) = truncated_threshold;
%     end
% end
%
% Normalize 0~1
%
% % % % % 
% for ii = 1:size(RD_map_noclutter,3) % 模擬次數
%     max_val=max(max(RD_map_noclutter(:,:,ii)));
%     min_val=min(min(RD_map_noclutter(:,:,ii)));
%     RD_map_normal(:,:,ii)=(RD_map_noclutter(:,:,ii)-min_val)/(max_val-min_val);
% end   
% % % % % 

%
% Dynamic Range Compression (DRC)
%
for ii = 1:size(RD_map_noclutter,3) % 模擬次數
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
% figure(1) % 只有target
% see1 = RD_map_label(:,:,1);
% mesh(see1,'edgecolor','r');
% 
% figure(2) % 只有noise
% see2 = reshape(RDmap_label_raw(1,:),N,M);
% mesh(see2,'edgecolor','r');
% zlim([0,10])
% 
% figure(3) % 沒有clutter
% see3 = reshape(RDmap_input_raw_noclutter2(1,:),N,M);
% mesh(see3,'edgecolor','r');
% 
% figure(10)
% imagesc(see3);
% 
% figure(4) % 經過truncated
% see4 = RD_map_truncated(:,:,1);
% mesh(see4,'edgecolor','r');
% title('Truncated')
% 
% figure(5) % Dynamic range compression
see5 = reshape(RDmap_input_raw_softknee(1, :), N, M);
% mesh(see5,'edgecolor','r');
% title('Dynamic range compression');
figure(11)
imagesc(see5);
% 
% figure(6) % Dynamic range compression
% see6 = reshape(RDmap_input_raw_log(1,:),N,M);
% mesh(see6,'edgecolor','r');
% title('Dynamic range compression: log')
% 
% figure(7) %  Normalize 0~1
% see7 = RD_map_normal(:,:,1);
% mesh(see7,'edgecolor','r');
% title('Normalize')

%
% validation
%
% save(['validation_noT_noclutter_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_noclutter'); % for CFAR
% save(['D:\YU\my_yolo3\DLCFAR\validation_noT_truncated_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_truncated'); % for DL-CFAR input
% save(['D:\YU\my_yolo3\DLCFAR\validation_noT_label_noise_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_label_raw'); % for DL-CFAR label(target)

% for jj=1:size(RD_map_softknee,3)
%     fileName = ['D:\YU\my_yolo3\VOC2007\JPEGImages_val\validation_noT_softknee_H' num2str(H) '_SNR' num2str(SNR) '_f' num2str(jj) '.mat'];
%     RD_map_yoloinput=RD_map_softknee(:,:,jj);
%     save(fileName,'RD_map_yoloinput'); % for YOLO-CFAR input
% end
% save(['validation_noT_label_target_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_label'); % for YOLO-CFAR label(target)

%
% training
%   create 5 folders CFAR_data, DL_input, DL_label, VOC2007/JPEGImages, YOLO_label at current directory
%   "CFAR_data" stores inputs for conventional CFAR detectors
%   "DL_input", "DL_label" stores training inputs and labels for DL_CFAR
%   "VOC2007/JPEGImages", "YOLO_label" stores training inputs and labels for YOLO_CFAR
%   notes that YOLO_input should store in \VOC2007\JPEGImages for convenient
% 
% save(['.\CFAR_data\CFAR_train_noT_noclutter_H',num2str(H),'_SNR',num2str(SNR),'.mat'],'RD_map_noclutter'); % for CFAR
% save(['.\DL_input\DL_input_train_noT_truncated_H',num2str(H),'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_truncated'); % for DL-CFAR input
% save(['.\DL_label\DL_label_train_noT_label_noise_H',num2str(H),'_SNR',num2str(SNR),'.mat'],'RDmap_label_raw'); % for DL-CFAR label(target)
% 
% for jj = 1:size(RD_map_softknee,3)
%     % % fileName = ['D:\BeginnerMatlabProjects\RDMap\VOC2007\JPEGImages\training_noT_softknee_H',num2str(H),'_SNR',num2str(SNR),'_f',num2str(jj),'.mat'];
%     fileName = ['.\VOC2007\JPEGImages\YOLO_input_train_noT_softknee_H',num2str(H),'_SNR',num2str(SNR),'_f',num2str(jj),'.mat'];
%     RD_map_yoloinput = RD_map_softknee(:,:,jj);
%     save(fileName,'RD_map_yoloinput'); % for YOLO-CFAR input
% end
% % % save(['training_noT_label_target_H',num2str(H),'_SNR',num2str(SNR),'.mat'],'RD_map_label'); % for YOLO-CFAR label(target)
% save(['.\YOLO_label\YOLO_label_train_noT_label_target_H',num2str(H),'_SNR',num2str(SNR),'.mat'],'RD_map_label');

%
% testing
%
% save(['testing_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_noclutter'); % for CFAR
% save(['D:\YU\my_yolo3\DLCFAR\testing_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_truncated'); % for DL-CFAR input
% save(['testing_label_target_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_label'); % label(target)
% SNR=20dB,single-target被蓋掉了?

% save(['testing_noT_noclutter_softVStrun_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_noclutter'); % for CFAR
% save(['D:\YU\my_yolo3\DLCFAR\testing_noT_soft_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_input_raw_softknee'); % for DL-CFAR input
% save(['D:\YU\my_yolo3\DLCFAR\testing_noT_label_noise_soft_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RDmap_label_raw'); % for DL-CFAR label(target)
% for jj=1:size(RD_map_truncated,3)
%     fileName = ['D:\YU\my_yolo3\mAP-master\input\images-optional(trun_20dB_H1)\f' num2str(jj) '_H' num2str(H) '_SNR' num2str(SNR) '.mat'];
%     RD_map_yoloinput=RD_map_truncated(:,:,jj);
%     save(fileName,'RD_map_yoloinput'); % for YOLO-CFAR input
% end
% save(['testing_noT_label_target_trun_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_label'); % for YOLO-CFAR label(target)
% save(['testing_trainDNN_trun_H' ,num2str(H) ,'_SNR',num2str(SNR),'.mat'],'RD_map_truncated'); % for DNN input

save para H SNR N M; % testing_trainDNN_softknee_20_H ???

% RD_map_yoloinput = see5;
% save(['D:\YU\my_yolo3\mAP-master\input\f1_H1_SNR6','.mat'],'RD_map_yoloinput');
% 
