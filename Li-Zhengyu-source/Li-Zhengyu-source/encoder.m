clear all;
close all;
clc;
I = imread('lena512.bmp');
I = im2double(I);
T = dctmtx(8); % dct matrix
%Performing DCT on blocks of 8 by 8
dct = @(block_struct) T * block_struct.data * T';
B = blockproc(I,[8 8],dct);
%B = ceil(B); 
% Quantization Matrix
% step=1;
Q=1;
 q_mtx =  Q* [16 11 10 16 24 40 51 61; 
            12 12 14 19 26 58 60 55;
            14 13 16 24 40 57 69 56; 
            14 17 22 29 51 87 80 62;
            18 22 37 56 68 109 103 77;
            24 35 55 64 81 104 113 92;
            49 64 78 87 103 121 120 101;
            72 92 95 98 112 100 103 99];
   %PErforming Quantization by Dividing with q_mtx on blocks of 8 by 8
   c = @(block_struct) (block_struct.data) ./ q_mtx;        
   B2 = blockproc(B,[8 8],c);
%      B2 = ceil(B2)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DCT Quantize Finished
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DC coefficients prediction Start
  temp = B2(1,1);
  for j = 1:8:512
     if j > 1
         temp0 = temp;
         temp = B2(1,j);
         B2(1,j)= B2(1,j)-temp0;
     end
  end
   temp0 = temp;
  for i=9:8:512
    for j = 1:8:512
          temp0 = temp;
          temp = B2(i,j);
          B2(i,j)= B2(i,j)-temp0;
     end
  end
% DC coefficients prediction finished
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Zig-zag scan start
 
% aa = reshape(B2,1,262144)
% zigzag = [ 0, 1, 8, 16, 9, 2, 3, 10, ...
%       17, 24, 32, 25, 18, 11, 4, 5, ...
%       12, 19, 26, 33, 40, 48, 41, 34, ...
%       27, 20, 13, 6, 7, 14, 21, 28, ...
%       35, 42, 49, 56, 57, 50, 43, 36, ...
%       29, 22, 15, 23, 30, 37, 44, 51, ...
%       58, 59, 52, 45, 38, 31, 39, 46, ...
%       53, 60, 61, 54, 47, 55, 62, 63];
% 
% zigzag = zigzag + 1;  % 
% aa = @(block_struct) reshape( (block_struct.data),1,64); % 
% B3 = blockproc(B2,[8 8],aa);
order = [
 0  1  5  6 14 15 27 28
 2  4  7 13 16 26 29 42
 3  8 12 17 25 30 41 43
 9 11 18 24 31 40 44 53
10 19 23 32 39 45 52 54
20 22 33 38 46 51 55 60
21 34 37 47 50 56 59 61
35 36 48 49 57 58 62 63
];
%aa=sortrows( [reshape(order, [], 1)]);
%aa =  @(block_struct)  sortrows( [reshape(order, [], 1)], 2);
aa =  @(block_struct)  sortrows( [reshape((block_struct.data), [], 1)  reshape(order, [], 1)], 2);
B3 = blockproc(B2,[8 8],aa);
%Zig zag scan end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 runl=0;
 y=1;
 z=1;
 s=1;
  for j=1:2:128
      for i=1:1:4096
             test(z,1) = B3(i,j);
            h(z,1)= ceil(test(z,1));
             if test(z,1) == 0
                 runl=runl+1;
                 z=z+1;
             else 
                 run(y,1)=runl;
                 if test(z,1) <= 1 && test(z,1) >= -1
                    s=1;
                elseif test(z,1)<=-2 &&test(z,1)>= -3 || test(z,1)>=2 && test(z,1) <= 3
                    s=2;
                elseif test(z,1)<=-4 &&test(z,1)>= -7 || test(z,1)>=4 && test(z,1) <= 7
                    s=3;
                elseif test(z,1)<=-8 &&test(z,1)>= -15 || test(z,1)>=8 && test(z,1) <= 15
                    s=4;
                elseif test(z,1)<=-16 && test(z,1)>=-31  || test(z,1)>=16 && test(z,1) <=31
                    s=5;
               elseif test(z,1)<=-32 && test(z,1)>=-63  || test(z,1)>=32 && test(z,1) <=63
                    s=6;
                    elseif test(z,1)<=-64 && test(z,1)>=-127  || test(z,1)>=64 && test(z,1) <=127
                        s=7;
                        elseif test(z,1)<=-128 && test(z,1)>=-255  || test(z,1)>=128 && test(z,1) <=255
                            s=8;
                            elseif test(z,1)<=-256 && test(z,1)>=-511  || test(z,1)>=256 && test(z,1) <=511
                                s=9;
                                elseif test(z,1)<=-512  || test(z,1)>=512 
                                    s=10;
                end
                 run(y,2)=s;
                 run(y,3)=h(z,1);
                 y=y+1;
                 z=z+1;
                 runl=0;
             end
          end
  end
       
  
   %%%% ru size,amp end; stored in run.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     x=[2 2 3 4 5];
%     tabulate(x(:))
 
%   [M,N] = size(run);
%   P = zeros(1,100);
% for i = 0:100
%      P(i+1) = length(find(run == i))/(M*N);
% end
% symbols = [1:256];
% dict = huffmandict(symbols,P);
% enco = huffmanenco(run,dict);
%  
% run3=[];
% for i=1:1:2
%     for j=1:1:3
% run3(i,j) = run(i,j);
%     end
% end
sum1=0;
[M,N] = size(run);  
symbols = -100:100;
for i = -100:100
      P(i+101) = length(find(abs(run) == i))/(M*N);
      sum1 = sum1 + P(i+101);  
end
dict = huffmandict(symbols,P);
 
run2 = reshape(run,1,M*N);
comp = huffmanenco(run2,dict);
 
%%%huffman code end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('C:\Users\zlr25\Documents\MATLAB\image.bit','wb') 
fwrite(fid,comp,'int')
fclose(fid)
 
 
 


