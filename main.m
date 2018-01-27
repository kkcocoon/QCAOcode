
%% matlab code for QCAO(Quadtree Coding with Adaptive scanning order and Optimized truncation)
% unoptimized, without head information, without entropy coding.
% 
% Reference:
% 
% Hui Liu, Ke-Kun Huang*, Chuan-Xian Ren, Yu-Feng Yu and Zhao-Rong Lai. 
% Quadtree Coding with Adaptive Scanning Order for Space-borne Image Compression,
% Signal Processing: Image Communication, Vol. 55, No. 7, pp. 1-9, 2017.
% 
% Email: kkcocoon@163.com
%

clc;clear;

%% Input 
imname = 'Detroit.bmp';
I_Orig = double(imread(imname));
[row, col] = size(I_Orig);
blksize = 64;  
n_log = log2(row); 
level = floor(n_log);
n_min = 1;
brates = [0.0313, 0.0625, 0.125, 0.25, 0.5, 1];  % The bit rates to calculate PSNR

I_Dec = wavecdf97(I_Orig, level);

%% Coding
[out_code, blklen, n_max, n_min, out_S,out_R,out_N] = encode(I_Dec, blksize, n_min);    

%% Decoding
disp([ 'aa_QTCAO_' imname(1:end-4) '=[']);
for rate=brates
    I_DecR = decode(out_code, blklen, n_max, n_min, blksize, row, rate, out_S,out_R,out_N);
    
    I_Rec = wavecdf97(I_DecR, -level);
    MSE = sum(sum((I_Rec - I_Orig).^2))/(row*row);
    PSNR = 10*log10(255*255/MSE);
    SSIM = fun_SSIM(I_Rec,I_Orig);
    disp([sprintf('%.4f',rate) ' ' sprintf('%.2f',PSNR) ' ' sprintf('%.4f',SSIM)]);   
end
disp('];');