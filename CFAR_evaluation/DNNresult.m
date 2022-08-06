clear all
% Soft
% close all
% % load testing_noT_clutter_CNR10_20_H1_SNR-4.mat
% load testing_noT_label_target_clutter_CNR10_20_H1_SNR-4.mat % 需改SNR

% % load testing_noT_noclutter_H4_SNR-4.mat
load testing_noT_label_target_trun_H1_SNR-4.mat % 需改SNR
load DNN_testing_input_label.mat
DNN_output_raw = csvread('testing_DNN_finalresult.csv');

for qqq = 1:size(DNN_output_raw,1)
    DNN_output(:,:,qqq) = reshape(DNN_output_raw(qqq,:),[3,3]);
end

filename = 'result.csv';
result = csvread(filename);
target_num = size(result,2)/5; % bounding box個數
for i=1:size(result,1)
    yolo_output(:,:,i)=reshape(result(i,:),[5,target_num])';
end
% % 1.找相鄰的框框(檢查第二三行是否一個一樣一個相差1) 2.confidence相加大約等於1(0.9~1.2)
% for iii=1:size(yolo_output,3)
%     v = (yolo_output(:,:,iii));
%         v(all(v == 0, 2),:) = [];
%     if size(v,1)~=1       
%     [c1,i1,j1]=unique(v(:,2));
%     if size(v,1)-length(c1)~=0 % 表示第二行有相同的值
%         
%         [M1,F1,C1] = mode(v(:,2)); % 找出第二行中出現次數最多的值
%         %ccc=zeros(,);
%         for cc = 1:length(C1{1}) % 須擴展到多目標
%             ccc1 = find(v(:,2)==C1{1}(cc)); % 得到第二行第幾列出現一樣的數字
%         end
%         if abs(v(ccc1(1),3)-v(ccc1(2),3))==1
%             if v(ccc1(1),1)+v(ccc1(2),1)>0.8 && v(ccc1(1),1)+v(ccc1(2),1)<1.4
%                 yolo_output(ccc1(1),1,iii)=1;
%                 yolo_output(ccc1(2),1,iii)=1;
%             end
%         end
%     end
%     [c2,i2,j2]=unique(v(:,3));
%     if size(v,1)-length(c2)~=0 % 表示第三行有相同的值
%         [M2,F2,C2] = mode(v(:,3)); % 找出第三行中出現次數最多的值
%         %ccc=zeros(,);
%         for cc = 1:length(C2{1}) % 須擴展到多目標
%             ccc2 = find(v(:,3)==C2{1}(cc)); % 得到第三行第幾列出現一樣的數字
%         end
%         if abs(v(ccc2(1),2)-v(ccc2(2),2))==1
%             if v(ccc2(1),1)+v(ccc2(2),1)>0.8 && v(ccc2(1),1)+v(ccc2(2),1)<1.4
%                 yolo_output(ccc2(1),1,iii)=1;
%                 yolo_output(ccc2(2),1,iii)=1;
%             end
%         end
%     end
%     else
%         yolo_output(1,1,iii)=1;
%     end
% end
for j=1:size(yolo_output,3) % 模擬次數
    RD_map_bboxnum(j,:)=sum(sum(yolo_output(:,1,j)~=0));
end

% iii為模擬次數
% ii為總共檢測出的bounding box個數
% jj為每張RDmap檢測出的bounding box
confidence_map=zeros(16,16,size(yolo_output,3));

iii=1;
ii=1;
while iii~=size(yolo_output,3)+1
    for jj=1:RD_map_bboxnum(iii)        
        if yolo_output(jj,3,iii)==0
            yolo_output(jj,3,iii)=yolo_output(jj,3,iii)+1;
            yolo_output(jj,5,iii)=yolo_output(jj,5,iii)+1;
        end
        if yolo_output(jj,2,iii)==0
            yolo_output(jj,2,iii)=yolo_output(jj,2,iii)+1;
            yolo_output(jj,4,iii)=yolo_output(jj,4,iii)+1;
        end
        if yolo_output(jj,3,iii)-yolo_output(jj,5,iii)==-1
            yolo_output(jj,5,iii)=yolo_output(jj,5,iii)+1;
        end
        if yolo_output(jj,2,iii)-yolo_output(jj,4,iii)==-1
            yolo_output(jj,4,iii)=yolo_output(jj,4,iii)+1;
        end
        if yolo_output(jj,3,iii)-yolo_output(jj,5,iii)==-3
            yolo_output(jj,5,iii)=yolo_output(jj,5,iii)-1;
        end
        if yolo_output(jj,2,iii)-yolo_output(jj,4,iii)==-3
            yolo_output(jj,4,iii)=yolo_output(jj,4,iii)-1;
        end
        v = (yolo_output(:,:,iii));
        v(all(v == 0, 2),:) = [];
        d = DNN_output(:,:,ii);
        if size(v,1) == 1
            d(d>0.9)=1;
%             d(d<0.5)=0;
            confidence_map(yolo_output(jj,3,iii):yolo_output(jj,5,iii),yolo_output(jj,2,iii):yolo_output(jj,4,iii),iii)=yolo_output(jj,1,iii)*d;
        else
            temp_map=zeros(16,16);
            temp_map(yolo_output(jj,3,iii):yolo_output(jj,5,iii),yolo_output(jj,2,iii):yolo_output(jj,4,iii))=yolo_output(jj,1,iii)*DNN_output(:,:,ii);
            a(:,:,jj)=DNN_output(:,:,ii);
            confidence_map(:,:,iii)=max(confidence_map(:,:,iii),temp_map);
        end
        ii=ii+1;
    end    
    iii=iii+1;
end
%%  畫出YOLO-CFAR的Pd Pfa圖
T_YOLO = (0.95:-0.05:0);% 0.05  0.95
for iiii=1:length(T_YOLO)
    TotalCount_D_YOLO = 0;
    SuccessCount_D_YOLO = 0;
    TotalCount_FA_YOLO = 0;
    FACount_FA_YOLO = 0;
    for TimeIndex = 1:size(RD_map_label,3) %36
        DetectionMap_YOLO = double(confidence_map(:,:,TimeIndex) >= T_YOLO(iiii));
        
%         testmap(:,:,TimeIndex)=DetectionMap_YOLO&RD_map_label(:,:,TimeIndex);
%         if sum(sum(testmap(:,:,TimeIndex))) ~=3
%             fprintf('%d,',TimeIndex)
%         end
        
        Detect_truthtable_YOLO = DetectionMap_YOLO(find(RD_map_label(:,:,TimeIndex) == 1));
        TotalCount_D_YOLO = TotalCount_D_YOLO + length(Detect_truthtable_YOLO);   
        SuccessCount_D_YOLO = SuccessCount_D_YOLO + length(find(Detect_truthtable_YOLO==1));
        
        [x,y] = find(RD_map_label(:,:,TimeIndex,1) == 1);
        
        for index = 1 : size(x,1)
%         if y(index) ~= 1 && x(index) ~= 1  
            DetectionMap_YOLO(x(index)-1,y(index)-1) = 0.5;
%         end
%         if x(index) ~= 1
            DetectionMap_YOLO(x(index)-1,y(index)) = 0.5;
%         end    
%         if x(index) ~= 1 && y(index) ~= 16
            DetectionMap_YOLO(x(index)-1,y(index)+1) = 0.5;
%         end
%         if y(index) ~= 1
            DetectionMap_YOLO(x(index),y(index)-1) = 0.5;
%         end            
%         if y(index) ~= 16
            DetectionMap_YOLO(x(index),y(index)+1) = 0.5;
%         end                
%         if x(index) ~= 16 && y(index) ~= 1
            DetectionMap_YOLO(x(index)+1,y(index)-1) = 0.5;
%         end                    
%         if x(index) ~= 16
            DetectionMap_YOLO(x(index)+1,y(index)) = 0.5;
%         end                        
%         if x(index) ~= 16 && y(index) ~= 16
            DetectionMap_YOLO(x(index)+1,y(index)+1) = 0.5;
%         end
        end
        
%         testmap(:,:,TimeIndex)=DetectionMap_YOLO-RD_map_label(:,:,TimeIndex);%false alarm
%         if length(find(testmap(:,:,TimeIndex)==1)) ~=0
%             fprintf('%d,',TimeIndex)
%         end
        
        FA_truthtable_YOLO = DetectionMap_YOLO(find(RD_map_label(:,:,TimeIndex) == 0));
        TotalCount_FA_YOLO = TotalCount_FA_YOLO + length(FA_truthtable_YOLO);
        FACount_FA_YOLO = FACount_FA_YOLO + length(find(FA_truthtable_YOLO==1));      
    end
    PDR_YOLO(iiii) = SuccessCount_D_YOLO/TotalCount_D_YOLO;
    FAR_YOLO(iiii) = FACount_FA_YOLO/TotalCount_FA_YOLO;
end

save yolo_roc PDR_YOLO FAR_YOLO
figure()
semilogx(FAR_YOLO,PDR_YOLO,'r-o')
xlabel('false alarm rate')
ylabel('detection rate')
ylim([0 1])
xlim([10^-8 10^0])
grid on
