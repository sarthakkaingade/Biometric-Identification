clc;
clear all;
close all;

[f,p] = uigetfile;
IMG = imread(cat(2,p,f));
[rows,columns] = size(IMG);
figure;
imshow(IMG);
a = mod(rows,10);
b = mod(columns,10);
if ( a ~= 0 )
    IMG(rows + (10 - a),:) = 0;
    rows = rows + (10 - a);
end
if ( b ~= 0 )
    IMG(:,columns + (10 - b)) = 0;
    columns = columns + (10 - b);
end
W = 10;
k = 0;
theta = zeros(rows / W, columns / W);
sobelx = [-1 -2 -1; 0 0 0; 1 2 1];
sobely = [-1 0 1; -2 0 2; -1 0 1];
for row = 1 : rows / W
    for col = 1 : columns / W
        img = IMG(((row-1) * W) + 1:((row-1) * W) + W,((col-1) * W) + 1:((col-1) * W) + W);
        Gx = filter2(sobelx,img);
        Gy = filter2(sobely,img);
        N = 0;
        for i = 2 : W -1
            for j = 2 : W - 1
                N = N + (2 * Gx(i,j) * Gy(i,j));
            end
        end
        D = 0;
        for i = 2 : W - 1
            for j = 2 : W - 1
                D = D + ((Gx(i,j) ^ 2) - (Gy(i,j) ^ 2));
            end
        end
        theta(row,col) = 0.5 * atan(N / D);
        if ( (theta(row,col) >= 0) && (N <= 0))
             k = 0;
        elseif ( (theta(row,col) < 0) && (N >= 0))
             k = 1;
        elseif ( ((theta(row,col) < 0) && (N < 0)) || ((theta(row,col) >= 0) && (N > 0)) )
             k = 0.5;        
        end
        theta(row,col) = theta(row,col) + (k * pi);
        thetafinal(((row-1) * W) + 1:((row-1) * W) + W,((col-1) * W) + 1:((col-1) * W) + W) = theta(row,col);
    end
end
drawOrientation(IMG, thetafinal);

rows = rows / W;
columns = columns / W;
im = double(zeros(rows + 2, columns + 2));
PI = zeros(rows + 2, columns + 2);
im(2 : rows + 1, 2 : columns + 1) = theta;
for i = 2:rows + 1
    for j = 2:columns + 1
        a = [im(i-1,j-1) im(i-1,j) im(i-1,j+1) im(i,j-1) im(i,j+1) im(i+1,j-1) im(i+1,j) im(i+1,j+1)];
%         a = sort(a);
        for k = 1:8
            if (k == 8)
                b = a(8) - a(1);
            else
                b = a(mod(k,8)) - a(k);
            end
            if (b > (pi/2))
                PI(i,j) = PI(i,j) + (b - pi);
            elseif (b < (-pi/2))
                PI(i,j) = PI(i,j) + (b + pi);
            else
                PI(i,j) = PI(i,j) + b;
            end
        end
        PI(i,j) = PI(i,j)/pi;
    end
end
