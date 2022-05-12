clear all;
% load testing_trainDNN_softknee_clutter_CNR10_20_H1_SNR-4.mat
% load testing_noT_label_target_clutter_CNR10_20_H1_SNR-4.mat % 惠эSNR

% load testing_trainDNN_trun_H1_SNR-4.mat
% load testing_noT_label_target_trun_H1_SNR-4.mat % 惠эSNR

load f1_H1_SNR6.mat

filename = 'result.csv';
result = csvread(filename);
target_num = size(result,2)/5; % bounding box计
for i=1:size(result,1)
    yolo_output(:,:,i)=reshape(result(i,:),[5,target_num])';
end
boxnum=0;
for j=1:size(yolo_output,3) % 家览Ω计
    num=sum(sum(yolo_output(:,1,j)~=0));
    boxnum=boxnum+num;
end
iiii=1;
for ii=1:size(yolo_output,3) % 家览Ω计
    yolo_temp=yolo_output([1:sum(sum(yolo_output(:,1,ii)~=0))],:,ii);    
    
    for iii=1:size(yolo_temp,1) % 块Bounding box计
        if yolo_temp(iii,3)==0
            yolo_temp(iii,3)=yolo_temp(iii,3)+1;
            yolo_temp(iii,5)=yolo_temp(iii,5)+1;
        end
        if yolo_temp(iii,2)==0
            yolo_temp(iii,2)=yolo_temp(iii,2)+1;
            yolo_temp(iii,4)=yolo_temp(iii,4)+1;
        end
        if yolo_temp(iii,3)-yolo_temp(iii,5)==-1
            yolo_temp(iii,5)=yolo_temp(iii,5)+1;
        end
        if yolo_temp(iii,2)-yolo_temp(iii,4)==-1
            yolo_temp(iii,4)=yolo_temp(iii,4)+1;
        end
        if yolo_temp(iii,3)-yolo_temp(iii,5)==-3
            yolo_temp(iii,5)=yolo_temp(iii,5)-1;
        end
        if yolo_temp(iii,2)-yolo_temp(iii,4)==-3
            yolo_temp(iii,4)=yolo_temp(iii,4)-1;
        end
        DNN_input(:,:,iiii)=RD_map_yoloinput([yolo_temp(iii,3):yolo_temp(iii,5)],[yolo_temp(iii,2):yolo_temp(iii,4)],ii);
%         DNN_input_label(:,:,iiii)=RD_map_label([yolo_temp(iii,3):yolo_temp(iii,5)],[yolo_temp(iii,2):yolo_temp(iii,4)],ii);
        iiii=iiii+1;
    end
end

for j = 1:size(DNN_input,3)% 家览Ω计
    DNN_input_raw(j,:) = reshape(DNN_input(:,:,j),[],9);
%     DNN_input_label_raw(j,:) = reshape(DNN_input_label(:,:,j),[],9);
end

% % train
% save(['D:\YU\my_yolo3\DNN\DNN_training_input.mat'],'DNN_input_raw');
% save(['D:\YU\my_yolo3\DNN\DNN_training_input_label.mat'],'DNN_input_label_raw');

% val
% save(['D:\YU\my_yolo3\DNN\DNN_val_input.mat'],'DNN_input_raw');
% save(['D:\YU\my_yolo3\DNN\DNN_val_input_label.mat'],'DNN_input_label_raw');

% test
save(['D:\YU\my_yolo3\DNN\DNN_testing_input_1.mat'],'DNN_input_raw');
% save(['DNN_testing_input_label.mat'],'DNN_input_label_raw');
% save(['D:\YU\my_yolo3\DNN\DNN_testing_input_label.mat'],'DNN_input_label_raw');
