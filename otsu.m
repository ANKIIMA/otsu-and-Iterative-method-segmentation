clear; clc;

%% OTSU
I=rgb2gray(imread('3.png'));
subplot(2, 2, 1)
imshow(I);
title('原图');

subplot(2,2,4);
imhist(I);
title('灰度直方图');

T = Otsu(double(I));     %使用大津法计算阈值


[m,n] = size(I);
for i = 1:m*n
    if I(i) <= T
        I(i) = 0;
    else
        I(i) = 255;
    end
end


%% 迭代法
subplot(2, 2, 2)
imshow(uint8(I));
title('OTSU');


A = rgb2gray(imread('3.png'));

t = repeat(A);

[m,n] = size(A);
for i = 1:m*n
    if A(i) < t
        A(i) = 0;
    else
        A(i) = 255;
    end
end
subplot(2,2,3);
imshow(uint8(A));
title('迭代法');



 

function value = Otsu(A)

T = min(A(:)):max(A(:));           
d = zeros(size(T));     
[m, n] = size(A);        
total = m*n;         

for i = 1 : length(T)

    iFg = 0;          % 前景
    iBg = 0;          % 背景
    FgSum = 0;    % 前景灰度和
    BgSum = 0;    % 背景灰度和
    for j = 1 : m
        for k = 1 : n
            temp = A(j, k);
            if temp > T(i)
                iFg = iFg + 1;   
                FgSum = FgSum + temp;
            else
                iBg = iBg + 1;
                BgSum = BgSum + temp;
            end
        end
    end
    wf = iFg/total;     % 前景比例
    wb = iBg/total;     % 背景比例
    uf = FgSum/iFg;     % 前景灰度平均值
    ub = BgSum/iBg;     % 背景灰度平均值
    d(i) = wf*wb*(uf - ub)*(uf - ub);     % 计算方差
end
[~, index] = max(d);             % 最大值下标
value = T(index);

end

function t = repeat(A)
t = mean2(double(A));

while true
    t1 = 0.5*(mean(A(A>t))+mean(A(A<t)));
    if abs(t - t1)<1e-11
        break
    end
    t = t1;
end

end

