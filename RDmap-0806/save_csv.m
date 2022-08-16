%
% generate train.csv and test.csv
%
fid1 = fopen(['D:\Datasets\RD_maps\', 'train.csv'], 'a');
fid2 = fopen(['D:\Datasets\RD_maps\', 'test.csv'], 'a');
max_iter = 2000; % total number
num_test = 1000;
num_train = max_iter - num_test;
for iter = 1:num_train
    % fprintf(fid1,'%d.png,%d.txt\n', iter, iter); % 16-by-16 gray scale images
    fprintf(fid1,'%d_sc.jpg,%d.txt\n', iter, iter); % 416-bi-416 scaled color images
end
fprintf('%d\n', num_train);
for iter = num_train + 1:max_iter
    % fprintf(fid2,'%d.png,%d.txt\n', iter, iter); % 16-by-16 gray scale images
    fprintf(fid2,'%d_sc.jpg,%d.txt\n', iter, iter); % 416-bi-416 scaled color images
end
fprintf('%d\n', max_iter);
fclose(fid1);
fclose(fid2);
