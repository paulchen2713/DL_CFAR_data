%
% generate train.csv, valid.csv and test.csv
%
fid1 = fopen(['D:\Datasets\RD_maps\', 'train.csv'], 'a');
fid2 = fopen(['D:\Datasets\RD_maps\', 'test.csv'], 'a');
max_iter = 1600; % total number
num_test = 100;
num_train = max_iter - num_test;
for iter = 1:num_train
    fprintf(fid1,'%d.jpg,%d.txt\n', iter, iter);
end
for iter = num_train + 1:max_iter
    fprintf(fid2,'%d.jpg,%d.txt\n', iter, iter);
end
fclose(fid1);
fclose(fid2);
