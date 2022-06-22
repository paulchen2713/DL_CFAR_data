%
% DNNinput main program
%
clear;
clc;
%
% 沒有load的資料 No such file or directory.
% 
% load testing_trainDNN_softknee_clutter_CNR10_20_H1_SNR-4.mat
% load testing_noT_label_target_clutter_CNR10_20_H1_SNR-4.mat % 需改SNR
%
% load testing_trainDNN_trun_H1_SNR-4.mat
% load testing_noT_label_target_trun_H1_SNR-4.mat % 需改SNR
%
% load f1_H1_SNR6.mat 
%
filename = 'result.csv';
result = csvread(filename);
target_num = size(result,2)/5; % bounding box 個數
for i = 1:size(result,1)
    yolo_output(:,:,i)=reshape(result(i,:),[5,target_num])';
end
boxnum = 0;
for j = 1:size(yolo_output,3) % 模擬次數
    num=sum(sum(yolo_output(:,1,j)~=0));
    boxnum=boxnum+num;
end
iiii=1;
for ii=1:size(yolo_output,3) % 模擬次數
    yolo_temp=yolo_output([1:sum(sum(yolo_output(:,1,ii)~=0))],:,ii);    
    %
    for iii=1:size(yolo_temp,1) % 輸出的Bounding box個數
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
        % DNN_input_label(:,:,iiii)=RD_map_label([yolo_temp(iii,3):yolo_temp(iii,5)],[yolo_temp(iii,2):yolo_temp(iii,4)],ii);
        iiii=iiii+1;
    end
end
%
for j = 1:size(DNN_input,3) % 模擬次數
    DNN_input_raw(j,:) = reshape(DNN_input(:,:,j),[],9);
    % DNN_input_label_raw(j,:) = reshape(DNN_input_label(:,:,j),[],9);
end
%
% train
%
% save(['D:\YU\my_yolo3\DNN\DNN_training_input.mat'],'DNN_input_raw');
% save(['D:\YU\my_yolo3\DNN\DNN_training_input_label.mat'],'DNN_input_label_raw');
%
% validation
%
% save(['D:\YU\my_yolo3\DNN\DNN_val_input.mat'],'DNN_input_raw');
% save(['D:\YU\my_yolo3\DNN\DNN_val_input_label.mat'],'DNN_input_label_raw');
%
% test
%
% save('D:\YU\my_yolo3\DNN\DNN_testing_input_1.mat','DNN_input_raw');
% save(['DNN_testing_input_label.mat'],'DNN_input_label_raw');
% save(['D:\YU\my_yolo3\DNN\DNN_testing_input_label.mat'],'DNN_input_label_raw');
% 

