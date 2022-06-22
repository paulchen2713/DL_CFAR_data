%
function [rdmap_clutter, rdmap_noclutter, label, noise] = newrdmap2(H, SNR)
    % fprintf('Y')
    %
    % Simulation settings
    %
    TotalSimulationTime = 1;    % number of RD maps
    CarrierFreq = 78*10^9;
    %
    % MIMO parameter settings
    %
    numTx = 1;
    numRx = 1;
    %
    % Frame-related parameters
    % RD map 大小為N*M
    N = 16;       % number of subcarrier
    M = 16;       % number of OFDM symbol
    Nfft  = N;    % number of FFT points in frequency domain
    frame = 1;    % number of frame used
    Mfft  = M;    % number of FFT points in time domain
    %
    % clutter parameters
    % 
    CNR = 15; %
    clutter_power = 0.5*10^(CNR/10);
    clutter_size_range = 5;
    clutter_size_Doppler = 5;
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
    d_unamb = c/2/SubcarrierSpacing;                  % unambiguity range
    v_unamb = c/2/CarrierFreq/PeriodOFDMsymbol_whole; % unambiguity velocity
    d_rel = d_unamb/N;    % search resolution of range
    v_rel = v_unamb/M;    % search resolution of velocity
    % 
    % Data matrix construction
    % 
    nmaxx = ones(H, 1);
    mmaxx = ones(H, 1);
    tn = zeros(H, 1);
    tm = zeros(H, 1);
    while sum(abs(nmaxx-tn)==1)>0 || sum(abs(mmaxx-tm)==1)>0
    %     fprintf('X')
    for g = 1: length(SNR)
        SNR_g = SNR(g);
        % target range setting
        RDmap_signal = cell(TotalSimulationTime,1); % 產生空細胞陣列
        RDmap_noise = cell(TotalSimulationTime,1);
        parameter_setting_record = [];
        for TimeIndex = 1:TotalSimulationTime
            tempMap = zeros(N,M);
            tempMap(1,1) = 1;
           while sum(tempMap(1,:))>0 || sum(tempMap(N,:))>0 || sum(tempMap(:,1))>0 || sum(tempMap(:,M))>0 || sum(sum(tempMap))<H
    %             fprintf('%d',sum(tempMap(1,:)~=0)>0 || sum(tempMap(N,:)~=0)>0 || sum(tempMap(:,1)~=0)>0 || sum(tempMap(:,M)~=0)>0);   
            for h = 1:H
                Range(h,1) = rand*d_unamb;
            end
            % target Doppler velocity setting
            for h = 1:H
                Vdop(h,1) = (2*rand-1)*v_unamb;
            end
    %         % target DoA setting
    %         for h = 1:H
    %             DoA(h,1) = (2*rand-1)*60; % FOV = 120 degrees
    %         end
            parameter_setting_record = [parameter_setting_record [Range;Vdop]]; 
            % transmitted signal (QPSK symbols)
            F_Tx = qammod(round(rand(N,M)*4+0.5,0)-1,4)/sqrt(2);  
            % channel effect
            F_Channel = cell(H,1);
            tempMap = zeros(N,M);
            target_Map(:,:,TimeIndex) = zeros(N,M);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            for index_target = 1:H
                NuiPhase = rand*2*pi;    % unpredictable phase difference between sources
                range=exp(-1j*2*pi*RoundTripTime(Range(index_target))*SubcarrierSpacing*[1:N].');
                doppler=exp(1j*2*pi*PeriodOFDMsymbol_whole*FreqDopp(Vdop(index_target))*[1:M]);
                F_Channel{index_target} =  F_Tx.*(range*doppler)*exp(1j*NuiPhase);
                RD_map_single_pure_target = abs(fft2(F_Channel{index_target}./F_Tx,N,M));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [nn,mm]=find(RD_map_single_pure_target == max(max(RD_map_single_pure_target)));
                tempMap(nn,mm) = 1; %%%%%%%%%%%%%%%%%%%%%
            end       
           end
            target_Map(:,:,TimeIndex) = tempMap;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % received signal
            F_Rx = cell(1,1);
            F_Rx_phase = cell(1,1);    
            F_Rx = zeros(N,M);        
            P_noise = 0.5;
            P_signal = P_noise*(10^(SNR_g/10));
            for index_target = 1:H
                F_Rx = F_Rx + sqrt(P_signal/2)*F_Channel{index_target};
            end
            % add sptially white noise
            P_Rx = mean(mean(F_Rx.*conj(F_Rx)));
            Z = sqrt(P_noise/2)*(randn(N,M)+1j*randn(N,M));
            F_Rx_n = F_Rx + Z;
            F_Rx_phase = F_Rx_n./F_Tx;
            Z_processed = Z./F_Tx;
            F_Rx_phase_signal_only = F_Rx./F_Tx;
            RDmap_signal{TimeIndex} = fft2(F_Rx_phase_signal_only,N,M)/sqrt(N)/sqrt(M);
            RDmap_noise{TimeIndex} = fft2(Z_processed,N,M)/sqrt(N)/sqrt(M);

            clutter_range_start = randi(N-clutter_size_range+1);
            clutter_range_end = clutter_range_start+clutter_size_range-1;
            clutter_Doppler_start = randi(M-clutter_size_Doppler+1);
            clutter_Doppler_end = clutter_Doppler_start+clutter_size_Doppler-1;
            clutter = zeros(N,M);
            clutter(clutter_range_start:clutter_range_end,clutter_Doppler_start:clutter_Doppler_end) = sqrt(clutter_power/2)+1j*sqrt(clutter_power/2);

            RDmap_full{TimeIndex} = RDmap_signal{TimeIndex} + RDmap_noise{TimeIndex}+clutter;
            RDmap_full_noclutter{TimeIndex} = RDmap_signal{TimeIndex} + RDmap_noise{TimeIndex};
        end
        % vectorization
        for TimeIndex = 1:TotalSimulationTime
            RDmap_input(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_full{TimeIndex}),1,N*M);
            RDmap_input_noclutter(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_full_noclutter{TimeIndex}),1,N*M);
            RDmap_label(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(abs(RDmap_noise{TimeIndex}),1,N*M); % noise
            RDmap_label_true_target(TimeIndex+(g-1)*TotalSimulationTime,1:N*M) = reshape(target_Map(:,:,TimeIndex),1,N*M); % target
        end
    end
    for i = 1 :  length(SNR)*TotalSimulationTime
        RDmap_input_raw(i,:) = RDmap_input(i,:).^2;
        RDmap_input_raw_noclutter(i,:) = RDmap_input_noclutter(i,:).^2;
        RDmap_label_raw(i,:) = RDmap_label(i,:).^2;
    end
    %% reshape N*M
    for ii = 1:size(RDmap_label_true_target,1)
        label(:,:,ii) = reshape(RDmap_label_true_target(ii,:),[N,M]); % 只有target
    end
    for ii = 1:size(RDmap_input_raw_noclutter,1)% size(A,1):A有幾列 = 模擬次數
        rdmap_clutter(:,:,ii) = reshape(RDmap_input_raw(ii,:),[N,M]); % target+noise+clutter
    end
    for ii = 1:size(RDmap_input_raw_noclutter,1)% size(A,1):A有幾列 = 模擬次數
        rdmap_noclutter(:,:,ii) = reshape(RDmap_input_raw_noclutter(ii,:),[N,M]); % target+noise
    end
    for ii = 1:size(RDmap_label_true_target,1)
        noise(:,:,ii) = RDmap_label_raw(ii,:); % 只有noise
    end
    %%
    % for ii = 1:size(RDmap_input_raw_noclutter,1)
    %     for h=1:H
    %     [nmax(h,ii),mmax(h,ii)] = find(rdmap(:,:,ii)==max(max(rdmap(:,:,ii))));
    %     [tn(h,ii),tm(h,ii)] = find(label(:,:,ii)==1);
    %     if nmax(h,ii)~=tn(h,ii) || mmax(h,ii)~=tm(h,ii)
    %         break;
    %     end
    %     end
    % end
    nmax=zeros(H,1);
    mmax=zeros(H,1);
    RD_map_noclutter_2=rdmap_noclutter;
    RD_map_noclutter_max=zeros(N,M);
    for h = 1:H
        [nmax(h),mmax(h)] = find(RD_map_noclutter_2(:,:)==max(max(RD_map_noclutter_2(:,:))));
        RD_map_noclutter_2(nmax(h),mmax(h)) = 0;
        RD_map_noclutter_max(nmax(h),mmax(h))=1;
    end   
    [nmaxx,mmaxx] = find(RD_map_noclutter_max(:,:)==1);
    [tn,tm] = find(label(:,:)==1);
    % if sum(abs(nmaxx-tn)==1)>0 || sum(abs(mmaxx-tm)==1)>0
    %    break;
    end   
end
