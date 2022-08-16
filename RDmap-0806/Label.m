% 
% Label main program
%   patch. June 22, 2022
% 
clear;
clc;
% 
% load training_noT_label_target_trun_clutter15_H4_SNR6.mat   % 需改SNR
% load validation_noT_label_target_trun_clutter15_H4_SNR6.mat % 需改SNR
% load training_noT_label_target_H4_SNR6.mat                  % 需改SNR
% load validation_noT_label_target_H4_SNR6.mat                % 需改SNR
% 
SNR = 6;
CNR = 15;
H = 4;
%
% para.mat 檔案 存有 H SNR N M 的資料 從 RD_map_dataset.m / RD_map_dataset_clutter.m 產生
% 但會覆寫掉上面 H SNR 設定 且 RD_map_dataset.m / RD_map_dataset_clutter.m 的 SNR 設定是像
% 0:5:10 跟 SNR = 某個定值 如6 衝突? 
load para.mat 
%
% fid=fopen(['D:\YU\my_yolo3\','2007_train.txt'],'a'); % 寫入檔案路徑(用"a"時，如果文字中已經存在資料，不會清空資料，而是在資料之後寫入，而"w"會清空原本的資料，重新寫入)
% fid=fopen(['D:\YU\my_yolo3\','2007_val.txt'],'a');
%
% 沒有 RD_map_label 的資料
% size(RD_map_label) = [N M length(SNR], e.g. 
for i = 1:size(RD_map_label, 3)
    [nn, mm] = find(RD_map_label(:, :, i) == 1); % target的位置
    for ii = 1:H
        ymin(ii) = nn(ii)-1; xmin(ii) = mm(ii)-1; % 左上
        ymax(ii) = nn(ii)+1; xmax(ii) = mm(ii)+1; % 右下
    end    
%     fprintf(fid,'VOC2007\\JPEGImages\\training_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1));
%     fprintf(fid,'VOC2007\\JPEGImages_val\\validation_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1));
% 
%     fprintf(fid,'VOC2007\\JPEGImages\\training_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2));
%     fprintf(fid,'VOC2007\\JPEGImages_val\\validation_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2));
% 
%     fprintf(fid,'VOC2007\\JPEGImages\\training_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3));
%     fprintf(fid,'VOC2007\\JPEGImages_val\\validation_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3));
% 
%     fprintf(fid,'VOC2007\\JPEGImages\\training_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3),xmin(4),ymin(4),xmax(4),ymax(4));
%     fprintf(fid,'VOC2007\\JPEGImages_val\\validation_noT_softknee_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3),xmin(4),ymin(4),xmax(4),ymax(4));
    %
    % clutter
    %
%     fprintf(fid,'VOC2007\\JPEGImages_trun\\training_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1));
%     fprintf(fid,'VOC2007\\JPEGImages_val_trun\\validation_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1));
% 
%     fprintf(fid,'VOC2007\\JPEGImages_trun\\training_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2));
%     fprintf(fid,'VOC2007\\JPEGImages_val_trun\\validation_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2));
% 
%     fprintf(fid,'VOC2007\\JPEGImages_trun\\training_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3));
%     fprintf(fid,'VOC2007\\JPEGImages_val_trun\\validation_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3));
% 
%     fprintf(fid,'VOC2007\\JPEGImages_trun\\training_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3),xmin(4),ymin(4),xmax(4),ymax(4));
%     fprintf(fid,'VOC2007\\JPEGImages_val_trun\\validation_noT_truncated_clutter%d_H%d_SNR%d_f%d.mat %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0 %d,%d,%d,%d,0\n',CNR,H,SNR,i,xmin(1),ymin(1),xmax(1),ymax(1),xmin(2),ymin(2),xmax(2),ymax(2),xmin(3),ymin(3),xmax(3),ymax(3),xmin(4),ymin(4),xmax(4),ymax(4));
end
% fclose(fid);

