clear; clc;

%% OTSU
I=rgb2gray(imread('3.png'));
subplot(2, 2, 1)
imshow(I);
title('ԭͼ');

subplot(2,2,4);
imhist(I);
title('�Ҷ�ֱ��ͼ');

T = Otsu(double(I));     %ʹ�ô�򷨼�����ֵ


[m,n] = size(I);
for i = 1:m*n
    if I(i) <= T
        I(i) = 0;
    else
        I(i) = 255;
    end
end


%% ������
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
title('������');



 

function value = Otsu(A)

T = min(A(:)):max(A(:));           
d = zeros(size(T));     
[m, n] = size(A);        
total = m*n;         

for i = 1 : length(T)

    iFg = 0;          % ǰ��
    iBg = 0;          % ����
    FgSum = 0;    % ǰ���ҶȺ�
    BgSum = 0;    % �����ҶȺ�
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
    wf = iFg/total;     % ǰ������
    wb = iBg/total;     % ��������
    uf = FgSum/iFg;     % ǰ���Ҷ�ƽ��ֵ
    ub = BgSum/iBg;     % �����Ҷ�ƽ��ֵ
    d(i) = wf*wb*(uf - ub)*(uf - ub);     % ���㷽��
end
[~, index] = max(d);             % ���ֵ�±�
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

