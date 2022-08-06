% 
% Label main program
%   patch. June 22, 2022
% 
clear;
clc;
% 
% load training_noT_label_target_trun_clutter15_H4_SNR6.mat   % �ݧ�SNR
% load validation_noT_label_target_trun_clutter15_H4_SNR6.mat % �ݧ�SNR
% load training_noT_label_target_H4_SNR6.mat                  % �ݧ�SNR
% load validation_noT_label_target_H4_SNR6.mat                % �ݧ�SNR
% 
SNR = 6;
CNR = 15;
H = 4;
%
% para.mat �ɮ� �s�� H SNR N M ����� �q RD_map_dataset.m / RD_map_dataset_clutter.m ����
% ���|�мg���W�� H SNR �]�w �B RD_map_dataset.m / RD_map_dataset_clutter.m �� SNR �]�w�O��
% 0:5:10 �� SNR = �Y�өw�� �p6 �Ĭ�? 
load para.mat 
%
% fid=fopen(['D:\YU\my_yolo3\','2007_train.txt'],'a'); % �g�J�ɮ׸��|(��"a"�ɡA�p�G��r���w�g�s�b��ơA���|�M�Ÿ�ơA�ӬO�b��Ƥ���g�J�A��"w"�|�M�ŭ쥻����ơA���s�g�J)
% fid=fopen(['D:\YU\my_yolo3\','2007_val.txt'],'a');
%
% �S�� RD_map_label �����
% size(RD_map_label) = [N M length(SNR], e.g. 
for i = 1:size(RD_map_label, 3)
    [nn, mm] = find(RD_map_label(:, :, i) == 1); % target����m
    for ii = 1:H
        ymin(ii) = nn(ii)-1; xmin(ii) = mm(ii)-1; % ���W
        ymax(ii) = nn(ii)+1; xmax(ii) = mm(ii)+1; % �k�U
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

