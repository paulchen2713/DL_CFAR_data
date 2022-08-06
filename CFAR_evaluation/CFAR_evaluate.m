clear all
% close all
% load testing_label_target_H3_SNR-4.mat % Τtarget
% load testing_H3_SNR-4.mat % ⊿Τ竒筁truncated⊿Τclutter

load testing_noT_label_target_20_H1_SNR-4.mat % Τtarget
load testing_noT_noclutter_20_H1_SNR-4.mat % ⊿Τ竒筁truncated⊿Τclutter
noise_Map_Deep_CFAR_unshaped = csvread('testing_DLCFAR_noise.csv');
% load para.mat
% load yolo_roc.mat 
N=16;M=16;

for ii = 1:size(noise_Map_Deep_CFAR_unshaped,1)
    noise_Map_Deep_CFAR(:,:,ii) = reshape(noise_Map_Deep_CFAR_unshaped(ii,:),[N,M]);
end

p = 0.1;
CA_detector = phased.CFARDetector2D('TrainingBandSize',[2,2], ...
    'ThresholdFactor','Auto','GuardBandSize',[1,1], ...
    'ProbabilityFalseAlarm',p,'Method','CA','NoisePowerOutputPort',true);
GOCA_detector = phased.CFARDetector2D('TrainingBandSize',[2,2], ...
    'ThresholdFactor','Auto','GuardBandSize',[1,1], ...
    'ProbabilityFalseAlarm',p,'Method','GOCA','NoisePowerOutputPort',true);
SOCA_detector = phased.CFARDetector2D('TrainingBandSize',[2,2], ...
    'ThresholdFactor','Auto','GuardBandSize',[1,1], ...
    'ProbabilityFalseAlarm',p,'Method','SOCA','NoisePowerOutputPort',true);
OS_detector_5 = phased.CFARDetector2D('TrainingBandSize',[2,2], ...
    'ThresholdFactor','Auto','GuardBandSize',[1,1], ...
    'ProbabilityFalseAlarm',p,'Method','OS','NoisePowerOutputPort',true,'Rank',5);
OS_detector_10 = phased.CFARDetector2D('TrainingBandSize',[2,2], ...
    'ThresholdFactor','Auto','GuardBandSize',[1,1], ...
    'ProbabilityFalseAlarm',p,'Method','OS','NoisePowerOutputPort',true,'Rank',10);
OS_detector_20 = phased.CFARDetector2D('TrainingBandSize',[2,2], ...
    'ThresholdFactor','Auto','GuardBandSize',[1,1], ...
    'ProbabilityFalseAlarm',p,'Method','OS','NoisePowerOutputPort',true,'Rank',20);
OS_detector_30 = phased.CFARDetector2D('TrainingBandSize',[2,2], ...
    'ThresholdFactor','Auto','GuardBandSize',[1,1], ...
    'ProbabilityFalseAlarm',p,'Method','OS','NoisePowerOutputPort',true,'Rank',30);
tic
% pad reverse part
L = 6; % L = window -1
Num_Padding = L/2;
for iii = 1:size(RD_map_noclutter,3) % 家览Ω计
    tempMatrix = [fliplr(RD_map_noclutter(:,end-Num_Padding+1:end,iii)) RD_map_noclutter(:,:,iii) fliplr(RD_map_noclutter(:,1:Num_Padding,iii))];
    tempMatrix = [flipud(tempMatrix(end-Num_Padding+1:end,:));tempMatrix;flipud(tempMatrix(1:Num_Padding,:))];
    RD_map_noclutter_pad(:,:,iii) = tempMatrix; % ﹍瓜(22*22)
end
for i = 1:size(noise_Map_Deep_CFAR,3)
    tempMatrix = [fliplr(noise_Map_Deep_CFAR(:,end-Num_Padding+1:end,i)) noise_Map_Deep_CFAR(:,:,i) fliplr(noise_Map_Deep_CFAR(:,1:Num_Padding,i))];
    tempMatrix = [flipud(tempMatrix(end-Num_Padding+1:end,:));tempMatrix;flipud(tempMatrix(1:Num_Padding,:))];
    noise_Map_Deep_CFAR_pad(:,:,i) = tempMatrix;
end

[N_ex] = size(RD_map_noclutter_pad(:,:,1),1); % 
[M_ex] = size(RD_map_noclutter_pad(:,:,1),2); % ︽
Ngc = CA_detector.GuardBandSize(2); % grardband︽计
Ngr = CA_detector.GuardBandSize(1); % grardband计
Ntc = CA_detector.TrainingBandSize(2);
Ntr = CA_detector.TrainingBandSize(1);
cutidx = []; % CUT
colstart = Ntc + Ngc + 1;
colend = N_ex - ( Ntc + Ngc);
rowstart = Ntr + Ngr + 1;
rowend = M_ex - ( Ntr + Ngr);
for m = colstart:colend
    for n = rowstart:rowend
        cutidx = [cutidx,[n;m]]; % CUT
    end
end
time_Deep = 0;
time_CA = 0;
time_GOCA = 0;
time_SOCA = 0;
time_OS_5 = 0;
time_OS_10 = 0;
time_OS_20 = 0;
time_OS_30 = 0;
% 璸衡noise level
for iii = 1:size(RD_map_noclutter_pad,3)
    tic
    [A temp] = CA_detector(noise_Map_Deep_CFAR_pad(:,:,iii),cutidx);
    mapCFAR_deep(:,:,iii) = reshape(temp,N,M);
    time_Deep = time_Deep + toc;
    tic
    [A temp] = CA_detector(RD_map_noclutter_pad(:,:,iii),cutidx);
    mapCFAR_CA(:,:,iii) = reshape(temp,N,M);
    time_CA = time_CA + toc;
    tic
    [A temp] = GOCA_detector(RD_map_noclutter_pad(:,:,iii),cutidx);
    mapCFAR_GOCA(:,:,iii) = reshape(temp,N,M);
    time_GOCA = time_GOCA + toc;
    tic
    [A temp] = SOCA_detector(RD_map_noclutter_pad(:,:,iii),cutidx);
    mapCFAR_SOCA(:,:,iii) = reshape(temp,N,M);
    time_SOCA = time_SOCA + toc;
    tic
    [A temp] = OS_detector_5(RD_map_noclutter_pad(:,:,iii),cutidx);
    mapCFAR_OS5(:,:,iii) = reshape(temp,N,M);
    time_OS_5 = time_OS_5 + toc;
    tic
    [A temp] = OS_detector_10(RD_map_noclutter_pad(:,:,iii),cutidx);
    mapCFAR_OS10(:,:,iii) = reshape(temp,N,M);
    time_OS_10 = time_OS_10 + toc;
    tic
    [A temp] = OS_detector_20(RD_map_noclutter_pad(:,:,iii),cutidx);
    mapCFAR_OS20(:,:,iii) = reshape(temp,N,M);
    time_OS_20 = time_OS_20 + toc;
    tic
    [A temp] = OS_detector_30(RD_map_noclutter_pad(:,:,iii),cutidx);
    mapCFAR_OS30(:,:,iii) = reshape(temp,N,M);
    time_OS_30 = time_OS_30 + toc;
end
time_Deep = time_Deep/size(RD_map_noclutter,3);
time_CA = time_CA/size(RD_map_noclutter,3);
time_GOCA = time_GOCA/size(RD_map_noclutter,3);
time_SOCA = time_SOCA/size(RD_map_noclutter,3);
time_OS_5 = time_OS_5/size(RD_map_noclutter,3);
time_OS_10 = time_OS_10/size(RD_map_noclutter,3);
time_OS_20 = time_OS_20/size(RD_map_noclutter,3);
time_OS_30 = time_OS_30/size(RD_map_noclutter,3);
toc
% 璸衡Pf Pd
P_FA_set = [1e-350 1e-300 1e-250 1e-200 1e-150 1e-100 1e-70 1e-60 1e-50 1e-40 1e-35 1e-30 1e-25 1e-22 1e-18 1e-14 1e-13 1e-10 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 0.5 0.9];
% P_FA_set = [1e-6];
for iiii = 1:length(P_FA_set)
    G =log(P_FA_set(iiii)^(-1)); %log(P_FA_set(iiii)^(-1)); P_FA_set(iiii)^(-1/40)-1
    TotalCount_D_DeepCFAR = 0;
    TotalCount_D_CACFAR = 0;
    TotalCount_D_GOCACFAR = 0;
    TotalCount_D_SOCACFAR = 0;
    TotalCount_D_OSCFAR_5 = 0;
    TotalCount_D_OSCFAR_10 = 0;
    TotalCount_D_OSCFAR_20 = 0;
    TotalCount_D_OSCFAR_30 = 0;    
    
    SuccessCount_D_DeepCFAR = 0;
    SuccessCount_D_CACFAR = 0;
    SuccessCount_D_GOCACFAR = 0;
    SuccessCount_D_SOCACFAR = 0;    
    SuccessCount_D_OSCFAR_5 = 0;
    SuccessCount_D_OSCFAR_10 = 0;
    SuccessCount_D_OSCFAR_20 = 0;
    SuccessCount_D_OSCFAR_30 = 0;
    
    TotalCount_FA_DeepCFAR = 0;
    TotalCount_FA_CACFAR = 0;
    TotalCount_FA_GOCACFAR = 0;
    TotalCount_FA_SOCACFAR = 0;
    TotalCount_FA_OSCFAR_5 = 0;
    TotalCount_FA_OSCFAR_10 = 0;
    TotalCount_FA_OSCFAR_20 = 0;
    TotalCount_FA_OSCFAR_30 = 0;    
    
    FACount_FA_DeepCFAR = 0;
    FACount_FA_CACFAR = 0;
    FACount_FA_GOCACFAR = 0;
    FACount_FA_SOCACFAR = 0;    
    FACount_FA_OSCFAR_5 = 0;
    FACount_FA_OSCFAR_10 = 0;
    FACount_FA_OSCFAR_20 = 0;
    FACount_FA_OSCFAR_30 = 0;
    
    for TimeIndex = 1:size(RD_map_noclutter,3)
        % threshold = G(scaling factor)*noiise level
        T_DeepCFAR = G*mapCFAR_deep(:,:,TimeIndex);
        T_CACFAR = G*mapCFAR_CA(:,:,TimeIndex);        
        T_GOCACFAR = G*mapCFAR_GOCA(:,:,TimeIndex);    
        T_SOCACFAR = G*mapCFAR_SOCA(:,:,TimeIndex);    
        T_OSCFAR_5 = G*(mapCFAR_OS5(:,:,TimeIndex)); 
        T_OSCFAR_10 = G*mapCFAR_OS10(:,:,TimeIndex);
        T_OSCFAR_20 = G*mapCFAR_OS20(:,:,TimeIndex);
        T_OSCFAR_30 = G*mapCFAR_OS30(:,:,TimeIndex);    
        %单T碞块1T块0
        DetectionMap_DLCFAR = double(RD_map_noclutter(:,:,TimeIndex) >= T_DeepCFAR);
        DetectionMap_CACFAR = double(RD_map_noclutter(:,:,TimeIndex) >= T_CACFAR);        
        DetectionMap_GOCACFAR = double(RD_map_noclutter(:,:,TimeIndex) >= T_GOCACFAR);        
        DetectionMap_SOCACFAR = double(RD_map_noclutter(:,:,TimeIndex) >= T_SOCACFAR);       
        DetectionMap_OSCFAR_5 = double(RD_map_noclutter(:,:,TimeIndex) >= T_OSCFAR_5);
        DetectionMap_OSCFAR_10 = double(RD_map_noclutter(:,:,TimeIndex) >= T_OSCFAR_10);
        DetectionMap_OSCFAR_20 = double(RD_map_noclutter(:,:,TimeIndex) >= T_OSCFAR_20);
        DetectionMap_OSCFAR_30 = double(RD_map_noclutter(:,:,TimeIndex) >= T_OSCFAR_30);    
        % detection rate pre-calculation
        % 絋粄targetΤ礚砆盎代 
        Detect_truthtable_DeepCFAR = DetectionMap_DLCFAR(find(RD_map_label(:,:,TimeIndex) == 1));
        Detect_truthtable_CACFAR = DetectionMap_CACFAR(find(RD_map_label(:,:,TimeIndex) == 1));        
        Detect_truthtable_GOCACFAR = DetectionMap_GOCACFAR(find(RD_map_label(:,:,TimeIndex) == 1));
        Detect_truthtable_SOCACFAR = DetectionMap_SOCACFAR(find(RD_map_label(:,:,TimeIndex) == 1));
        Detect_truthtable_OSCFAR_5 = DetectionMap_OSCFAR_5(find(RD_map_label(:,:,TimeIndex) == 1)); % RD_map_labelΤtarget
        Detect_truthtable_OSCFAR_10 = DetectionMap_OSCFAR_10(find(RD_map_label(:,:,TimeIndex) == 1));
        Detect_truthtable_OSCFAR_20 = DetectionMap_OSCFAR_20(find(RD_map_label(:,:,TimeIndex) == 1));
        Detect_truthtable_OSCFAR_30 = DetectionMap_OSCFAR_30(find(RD_map_label(:,:,TimeIndex) == 1));
        
        TotalCount_D_DeepCFAR = TotalCount_D_DeepCFAR + length(Detect_truthtable_DeepCFAR);
        TotalCount_D_CACFAR = TotalCount_D_CACFAR + length(Detect_truthtable_CACFAR);
        TotalCount_D_GOCACFAR = TotalCount_D_GOCACFAR + length(Detect_truthtable_GOCACFAR);
        TotalCount_D_SOCACFAR = TotalCount_D_SOCACFAR + length(Detect_truthtable_SOCACFAR);
        TotalCount_D_OSCFAR_5 = TotalCount_D_OSCFAR_5 + length(Detect_truthtable_OSCFAR_5);
        TotalCount_D_OSCFAR_10 = TotalCount_D_OSCFAR_10 + length(Detect_truthtable_OSCFAR_10);
        TotalCount_D_OSCFAR_20 = TotalCount_D_OSCFAR_20 + length(Detect_truthtable_OSCFAR_20);
        TotalCount_D_OSCFAR_30 = TotalCount_D_OSCFAR_30 + length(Detect_truthtable_OSCFAR_30);
        
        SuccessCount_D_DeepCFAR = SuccessCount_D_DeepCFAR + length(find(Detect_truthtable_DeepCFAR==1));
        SuccessCount_D_CACFAR = SuccessCount_D_CACFAR + length(find(Detect_truthtable_CACFAR==1));
        SuccessCount_D_GOCACFAR = SuccessCount_D_GOCACFAR + length(find(Detect_truthtable_GOCACFAR==1));
        SuccessCount_D_SOCACFAR = SuccessCount_D_SOCACFAR + length(find(Detect_truthtable_SOCACFAR==1));
        SuccessCount_D_OSCFAR_5 = SuccessCount_D_OSCFAR_5 + length(find(Detect_truthtable_OSCFAR_5==1));       
        SuccessCount_D_OSCFAR_10 = SuccessCount_D_OSCFAR_10 + length(find(Detect_truthtable_OSCFAR_10==1));
        SuccessCount_D_OSCFAR_20 = SuccessCount_D_OSCFAR_20 + length(find(Detect_truthtable_OSCFAR_20==1));
        SuccessCount_D_OSCFAR_30 = SuccessCount_D_OSCFAR_30 + length(find(Detect_truthtable_OSCFAR_30==1));
        % false alarm rate pre-calculation
        [x,y] = find(RD_map_label(:,:,TimeIndex,1) == 1);
        for index = 1 : size(x,1)
%         if y(index) ~= 1 && x(index) ~= 1  
            DetectionMap_DLCFAR(x(index)-1,y(index)-1) = 0.5;
            DetectionMap_CACFAR(x(index)-1,y(index)-1) = 0.5;
            DetectionMap_GOCACFAR(x(index)-1,y(index)-1) = 0.5;
            DetectionMap_SOCACFAR(x(index)-1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_5(x(index)-1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_10(x(index)-1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_20(x(index)-1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_30(x(index)-1,y(index)-1) = 0.5;            
%         end
%         if x(index) ~= 1
            DetectionMap_DLCFAR(x(index)-1,y(index)) = 0.5;
            DetectionMap_CACFAR(x(index)-1,y(index)) = 0.5;
            DetectionMap_GOCACFAR(x(index)-1,y(index)) = 0.5;
            DetectionMap_SOCACFAR(x(index)-1,y(index)) = 0.5;
            DetectionMap_OSCFAR_5(x(index)-1,y(index)) = 0.5;
            DetectionMap_OSCFAR_10(x(index)-1,y(index)) = 0.5;
            DetectionMap_OSCFAR_20(x(index)-1,y(index)) = 0.5;
            DetectionMap_OSCFAR_30(x(index)-1,y(index)) = 0.5;
%         end    
%         if x(index) ~= 1 && y(index) ~= 16
            DetectionMap_DLCFAR(x(index)-1,y(index)+1) = 0.5;
            DetectionMap_CACFAR(x(index)-1,y(index)+1) = 0.5;
            DetectionMap_GOCACFAR(x(index)-1,y(index)+1) = 0.5;
            DetectionMap_SOCACFAR(x(index)-1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_5(x(index)-1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_10(x(index)-1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_20(x(index)-1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_30(x(index)-1,y(index)+1) = 0.5;
%         end
%         if y(index) ~= 1
            DetectionMap_DLCFAR(x(index),y(index)-1) = 0.5;
            DetectionMap_CACFAR(x(index),y(index)-1) = 0.5;
            DetectionMap_GOCACFAR(x(index),y(index)-1) = 0.5;
            DetectionMap_SOCACFAR(x(index),y(index)-1) = 0.5;
            DetectionMap_OSCFAR_5(x(index),y(index)-1) = 0.5;
            DetectionMap_OSCFAR_10(x(index),y(index)-1) = 0.5;
            DetectionMap_OSCFAR_20(x(index),y(index)-1) = 0.5;
            DetectionMap_OSCFAR_30(x(index),y(index)-1) = 0.5;
%         end            
%         if y(index) ~= 16
            DetectionMap_DLCFAR(x(index),y(index)+1) = 0.5;
            DetectionMap_CACFAR(x(index),y(index)+1) = 0.5;
            DetectionMap_GOCACFAR(x(index),y(index)+1) = 0.5;
            DetectionMap_SOCACFAR(x(index),y(index)+1) = 0.5;
            DetectionMap_OSCFAR_5(x(index),y(index)+1) = 0.5;
            DetectionMap_OSCFAR_10(x(index),y(index)+1) = 0.5;
            DetectionMap_OSCFAR_20(x(index),y(index)+1) = 0.5;
            DetectionMap_OSCFAR_30(x(index),y(index)+1) = 0.5;
%         end                
%         if x(index) ~= 16 && y(index) ~= 1
            DetectionMap_DLCFAR(x(index)+1,y(index)-1) = 0.5;
            DetectionMap_CACFAR(x(index)+1,y(index)-1) = 0.5;
            DetectionMap_GOCACFAR(x(index)+1,y(index)-1) = 0.5;
            DetectionMap_SOCACFAR(x(index)+1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_5(x(index)+1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_10(x(index)+1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_20(x(index)+1,y(index)-1) = 0.5;
            DetectionMap_OSCFAR_30(x(index)+1,y(index)-1) = 0.5;
%         end                    
%         if x(index) ~= 16
            DetectionMap_DLCFAR(x(index)+1,y(index)) = 0.5;
            DetectionMap_CACFAR(x(index)+1,y(index)) = 0.5;
            DetectionMap_GOCACFAR(x(index)+1,y(index)) = 0.5;
            DetectionMap_SOCACFAR(x(index)+1,y(index)) = 0.5;
            DetectionMap_OSCFAR_5(x(index)+1,y(index)) = 0.5;
            DetectionMap_OSCFAR_10(x(index)+1,y(index)) = 0.5;
            DetectionMap_OSCFAR_20(x(index)+1,y(index)) = 0.5;
            DetectionMap_OSCFAR_30(x(index)+1,y(index)) = 0.5;
%         end                        
%         if x(index) ~= 16 && y(index) ~= 16
            DetectionMap_DLCFAR(x(index)+1,y(index)+1) = 0.5;
            DetectionMap_CACFAR(x(index)+1,y(index)+1) = 0.5;
            DetectionMap_GOCACFAR(x(index)+1,y(index)+1) = 0.5;
            DetectionMap_SOCACFAR(x(index)+1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_5(x(index)+1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_10(x(index)+1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_20(x(index)+1,y(index)+1) = 0.5;
            DetectionMap_OSCFAR_30(x(index)+1,y(index)+1) = 0.5;
%         end
        end
        % 单1よボΤfalse alarm
        FA_truthtable_DeepCFAR = DetectionMap_DLCFAR(find(RD_map_label(:,:,TimeIndex) == 0));
        FA_truthtable_CACFAR = DetectionMap_CACFAR(find(RD_map_label(:,:,TimeIndex) == 0));
        FA_truthtable_GOCACFAR = DetectionMap_GOCACFAR(find(RD_map_label(:,:,TimeIndex) == 0));
        FA_truthtable_SOCACFAR = DetectionMap_SOCACFAR(find(RD_map_label(:,:,TimeIndex) == 0));
        FA_truthtable_OSCFAR_5 = DetectionMap_OSCFAR_5(find(RD_map_label(:,:,TimeIndex) == 0));      
        FA_truthtable_OSCFAR_10 = DetectionMap_OSCFAR_10(find(RD_map_label(:,:,TimeIndex) == 0));
        FA_truthtable_OSCFAR_20 = DetectionMap_OSCFAR_20(find(RD_map_label(:,:,TimeIndex) == 0));
        FA_truthtable_OSCFAR_30 = DetectionMap_OSCFAR_30(find(RD_map_label(:,:,TimeIndex) == 0));

        TotalCount_FA_DeepCFAR = TotalCount_FA_DeepCFAR + length(FA_truthtable_DeepCFAR);
        TotalCount_FA_CACFAR = TotalCount_FA_CACFAR + length(FA_truthtable_CACFAR);
        TotalCount_FA_GOCACFAR = TotalCount_FA_GOCACFAR + length(FA_truthtable_GOCACFAR);
        TotalCount_FA_SOCACFAR = TotalCount_FA_SOCACFAR + length(FA_truthtable_SOCACFAR);
        TotalCount_FA_OSCFAR_5 = TotalCount_FA_OSCFAR_5 + length(FA_truthtable_OSCFAR_5);        
        TotalCount_FA_OSCFAR_10 = TotalCount_FA_OSCFAR_10 + length(FA_truthtable_OSCFAR_10);
        TotalCount_FA_OSCFAR_20 = TotalCount_FA_OSCFAR_20 + length(FA_truthtable_OSCFAR_20);
        TotalCount_FA_OSCFAR_30 = TotalCount_FA_OSCFAR_30 + length(FA_truthtable_OSCFAR_30);
        
        FACount_FA_DeepCFAR = FACount_FA_DeepCFAR + length(find(FA_truthtable_DeepCFAR==1));
        FACount_FA_CACFAR = FACount_FA_CACFAR + length(find(FA_truthtable_CACFAR==1));
        FACount_FA_GOCACFAR = FACount_FA_GOCACFAR + length(find(FA_truthtable_GOCACFAR==1));
        FACount_FA_SOCACFAR = FACount_FA_SOCACFAR + length(find(FA_truthtable_SOCACFAR==1));
        FACount_FA_OSCFAR_5 = FACount_FA_OSCFAR_5 + length(find(FA_truthtable_OSCFAR_5==1));
        FACount_FA_OSCFAR_10 = FACount_FA_OSCFAR_10 + length(find(FA_truthtable_OSCFAR_10==1));
        FACount_FA_OSCFAR_20 = FACount_FA_OSCFAR_20 + length(find(FA_truthtable_OSCFAR_20==1));
        FACount_FA_OSCFAR_30 = FACount_FA_OSCFAR_30 + length(find(FA_truthtable_OSCFAR_30==1));
    end
    % detection rate calculation
    PDR_DeepCFAR(iiii) = SuccessCount_D_DeepCFAR/TotalCount_D_DeepCFAR;
    PDR_CACFAR(iiii) = SuccessCount_D_CACFAR/TotalCount_D_CACFAR;
    PDR_GOCACFAR(iiii) = SuccessCount_D_GOCACFAR/TotalCount_D_GOCACFAR;
    PDR_SOCACFAR(iiii) = SuccessCount_D_SOCACFAR/TotalCount_D_SOCACFAR;
    PDR_OSCFAR_5(iiii) = SuccessCount_D_OSCFAR_5/TotalCount_D_OSCFAR_5;
    PDR_OSCFAR_10(iiii) = SuccessCount_D_OSCFAR_10/TotalCount_D_OSCFAR_10;
    PDR_OSCFAR_20(iiii) = SuccessCount_D_OSCFAR_20/TotalCount_D_OSCFAR_20;
    PDR_OSCFAR_30(iiii) = SuccessCount_D_OSCFAR_30/TotalCount_D_OSCFAR_30;
    % false alarm rate calculation
    FAR_DeepCFAR(iiii) = FACount_FA_DeepCFAR/TotalCount_FA_DeepCFAR;
    FAR_CACFAR(iiii) = FACount_FA_CACFAR/TotalCount_FA_CACFAR;
    FAR_GOCACFAR(iiii) = FACount_FA_GOCACFAR/TotalCount_FA_GOCACFAR;
    FAR_SOCACFAR(iiii) = FACount_FA_SOCACFAR/TotalCount_FA_SOCACFAR;
    FAR_OSCFAR_5(iiii) = FACount_FA_OSCFAR_5/TotalCount_FA_OSCFAR_5;
    FAR_OSCFAR_10(iiii) = FACount_FA_OSCFAR_10/TotalCount_FA_OSCFAR_10;
    FAR_OSCFAR_20(iiii) = FACount_FA_OSCFAR_20/TotalCount_FA_OSCFAR_20;
    FAR_OSCFAR_30(iiii) = FACount_FA_OSCFAR_30/TotalCount_FA_OSCFAR_30;
end
% save SNR20dB_H2_cfar_roc PDR_DeepCFAR FAR_DeepCFAR PDR_CACFAR PDR_GOCACFAR PDR_SOCACFAR PDR_OSCFAR_5 PDR_OSCFAR_10 PDR_OSCFAR_20 PDR_OSCFAR_30 FAR_CACFAR FAR_GOCACFAR FAR_SOCACFAR FAR_OSCFAR_5 FAR_OSCFAR_10 FAR_OSCFAR_20 FAR_OSCFAR_30
% open('SNR20dB_H2_w6(clutter)_DNNw5_200000_soft.fig');
% h_line=get(gca,'Children');%get linehandles
% xdata=get(h_line,'Xdata');!
% ydata=get(h_line,'Ydata');
figure()
% semilogx(xdata,ydata,'r-p')
% hold on
semilogx(FAR_DeepCFAR,PDR_DeepCFAR,'b-<')
hold on
semilogx(FAR_CACFAR,PDR_CACFAR,'m-o')
hold on
semilogx(FAR_GOCACFAR,PDR_GOCACFAR,'c->')
hold on
semilogx(FAR_SOCACFAR,PDR_SOCACFAR,'g-s')
hold on
semilogx(FAR_OSCFAR_5,PDR_OSCFAR_5,'k-h')
hold on
semilogx(FAR_OSCFAR_10,PDR_OSCFAR_10,'k-*')
hold on
semilogx(FAR_OSCFAR_20,PDR_OSCFAR_20,'k-^')
hold on
semilogx(FAR_OSCFAR_30,PDR_OSCFAR_30,'k-')
xlabel('false alarm rate')
ylabel('detection rate')
legend('DL-CFAR','CA-CFAR','GOCA-CFAR','SOCA-CFAR','OS-CFAR-5','OS-CFAR-10','OS-CFAR-20','OS-CFAR-30')
%'YOLO-CFAR',
grid on
ylim([0 1])
axis square     

% figure() % ⊿Τclutter
% mesh(RD_map_noclutter,'edgecolor','k');
% hold on
% mesh(T_DeepCFAR,'edgecolor','b');    
% 
% figure()
% mesh(RD_map_noclutter,'edgecolor','k');
% hold on
% mesh(T_CACFAR,'edgecolor','m');
% 
% figure()
% mesh(RD_map_noclutter,'edgecolor','k');
% hold on
% mesh(T_GOCACFAR,'edgecolor','c');  
% 
% figure()
% mesh(RD_map_noclutter,'edgecolor','k');
% hold on
% mesh(T_SOCACFAR,'edgecolor','g');
% 
% figure()
% mesh(RD_map_noclutter,'edgecolor','k');
% hold on
% mesh(T_OSCFAR_30,'edgecolor',[0.8 0.5 0]);  
        
        
        
 