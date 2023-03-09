% 
% convert .mat matrices into .jpg images
% 
% clear;
clc;

% NOTE. for "how to import or export a sequence of files" please refer to 
% https://www.mathworks.com/help/matlab/import_export/process-a-sequence-of-files.html

% Script seetings
load_data = 1;      % 1 means "do load"; 0 means "don't load"
save_imagesc = 1;   % 1 means "do save scaled color image"; 0 means "don't"
save_imwrite = 0;   % 1 means "do save gray scale image"; 0 means "don't"
total_num_file = 7193; % there are in total 7193 .mat files / .jpg images

% We store all the RD matrices in a 1-D cell array, which should end up 
% being a 1-by-7193 array. Each element within is a 256-by-64 RD matrix.
buffer = cell(1, total_num_file); 

% Loading all the .mat files might take a while, so if we load them once 
% and keep them in the workspace, we should be good.
if load_data == 1
    for i = 1:total_num_file
        fprintf('%d\n', i); % print out the current progress
        % set the .mat file path 
        file_path = ['D:/Datasets/RADA/RD_JPG/mats/',num2str(i),'.mat'];
        buffer{i} = importdata(file_path); 
    end
end

% Plotting and saving images in scaled color and grayscale.
for i = 1:total_num_file
    fprintf('%d\n', i); % print out the current progress
    
    % plot the image in scaled color mode
    figure = imagesc(buffer{i});
    set(gca,'XTick',[]) % remove the ticks in the x axis
    set(gca,'YTick',[]) % remove the ticks in the y axis
    set(gca,'Position', [0 0 1 1]) % make the axes occupy the hole figure
    
    % store the scaled color image as .jpg in default size, not 256-by-64
    curr_image = ['D:/Datasets/RADA/RD_JPG/imagesc/',num2str(i),'.jpg'];
    if save_imagesc == 1
        % fprintf('%s\n', curr_image);
        saveas(gcf, curr_image, 'jpg'); 
    end
    if save_imwrite == 1
        % store the grayscale image as .jpg in 256-by-64
        temp = rescale(buffer{i}, 0, 255); 
        temp = uint8(temp);
        curr_image = ['D:/Datasets/RADA/RD_JPG/imwrite/',num2str(i),'.jpg'];
        imwrite(temp, curr_image);
    end
end


