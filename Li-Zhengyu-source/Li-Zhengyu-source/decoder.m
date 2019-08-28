dsig = huffmandeco(comp,dict);
 a=0;b=0;
  cout=[];
 
for L=1:2:128   
  for k=1:64:4096
    j=1;
    in=[];
   for i=0:1:63
     in(j,1) = B3((i+k),L);
     j=j+1;
   end
tot_elem=length(in);
if tot_elem~=8*8
    error('Matrix dimensions do not coincide');
end 
% Initialise the output matrix
out=zeros(8,8);
 
cur_row=1;  cur_col=1;  cur_index=1;
while cur_index<=tot_elem
    
    if cur_row==1 & mod(cur_row+cur_col,2)==0 & cur_col~=8
        out(cur_row,cur_col)=in(cur_index);
        cur_col=cur_col+1;                          %move right at the top
        cur_index=cur_index+1;
        
    elseif cur_row==8 & mod(cur_row+cur_col,2)~=0 & cur_col~=8
        out(cur_row,cur_col)=in(cur_index);
        cur_col=cur_col+1;                          %move right at the bottom
        cur_index=cur_index+1;
        
    elseif cur_col==1 & mod(cur_row+cur_col,2)~=0 & cur_row~=8
        out(cur_row,cur_col)=in(cur_index);
        cur_row=cur_row+1;                          %move down at the left
        cur_index=cur_index+1;
        
    elseif cur_col==8 & mod(cur_row+cur_col,2)==0 & cur_row~=8
        out(cur_row,cur_col)=in(cur_index);
        cur_row=cur_row+1;                          %move down at the right
        cur_index=cur_index+1;
        
    elseif cur_col~=1 & cur_row~=8 & mod(cur_row+cur_col,2)~=0
        out(cur_row,cur_col)=in(cur_index);
        cur_row=cur_row+1;      cur_col=cur_col-1;  %move diagonally left down
        cur_index=cur_index+1;
        
    elseif cur_row~=1 & cur_col~=8 & mod(cur_row+cur_col,2)==0
        out(cur_row,cur_col)=in(cur_index);
        cur_row=cur_row-1;      cur_col=cur_col+1;  %move diagonally right up
        cur_index=cur_index+1;
        
    elseif cur_index==tot_elem                      %input the bottom right element
        out(end)=in(end);                           %end of the operation
        break                                       %terminate the operation
    end  
end
for i=1:1:8          
    for j=1:1:8      
cout(a+i,j+b)=  out(i,j);
    end
end
a=a+8;
if a == 512
    b=b+8;
    a=0;
end
end
end
%Inv zig zag end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = cout(1,1)
  for j = 9:8:512
         cout(1,j)= cout(1,j)+temp;
%          temp0 = temp;
         temp = cout(1,j);     
  end
  for i=9:8:512
    for j = 1:8:512
          cout(i,j) = cout(i,j)+temp;
          temp = cout(i,j);
     end
  end
  %%recontruct DC coefficients end
  %%%%%%%%%%%%%%%%%%%
Q=1.3;
 q_mtx =  Q* [16 11 10 16 24 40 51 61; 
            12 12 14 19 26 58 60 55;
            14 13 16 24 40 57 69 56; 
            14 17 22 29 51 87 80 62;
            18 22 37 56 68 109 103 77;
            24 35 55 64 81 104 113 92;
            49 64 78 87 103 121 120 101;
            72 92 95 98 112 100 103 99];
   c = @(block_struct) (block_struct.data) .* q_mtx;        
   B5 = blockproc(cout,[8 8],c);
invdct = @(block_struct) T' * block_struct.data * T;
I2 = blockproc(B5,[8 8],invdct);
imshow(I)
figure
imshow(I2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MSE
MSE=0;
for i=1:1:512
    for j=1:1:512
        D1 = (abs(I(i,j)-I2(i,j)));
        D = D1^2/512/512;
        MSE = MSE + D;   
    end
end
MSE = MSE * 100;
PSNR = 10*log10((255*255)/MSE);
CR = length(comp)/length(run);
figure
x=-1.7779:3.1647;
y=-3.253*x+49;
plot(x,y)
figure
t=3.1647:6.4126;
Z=-6.0206*t+35;
plot(t,Z)

