% this is the code for the implementation Histogram modification

%step1: clear previous data
clc;
clear all;
close all;
addpath('support')

%Step2: select an image for contrast enhancemnet
I=uigetfile('*.*','select the source image');
I=imread(I);
I=imresize(I,[512 512]);
S=I;
[M1 N1 K]=size(I);
figure,imshow(I);title('Original Image');
figure,imhist(I(:,:,1));title(['Histogram of Input']);
%-------------------------------------------------------
% step3: Apply Existing Histogram Equalization

outHE1=histeq(I(:,:,1));
outHE2=histeq(I(:,:,2));
outHE3=histeq(I(:,:,3));
outHE=cat(3,outHE1,outHE2,outHE3)/3;
figure,imshow(outHE),title('Histogram Equalization')
figure,imhist(outHE(:,:,1));title(['Histogram of HE']);
%-------------------------------------------------------
%step4: Proposed Method Modified HE 
b=20;
w=10;
alpha=2;
f=I;
% step5: Initialize variables
k = 0;
count = 0;
h = zeros(1, 256);
Threshold=128;

% Get the size of the image
[M, N] = size(f);

%step6: Loop through each pixel in the image
% for modification
for m = 1:M
    for n = 1:N
        % Calculate k
        k = k + mod(f(m, n) - f(m, max(1, n - 2)), 256);
        
        % Check condition and update histogram
        if k + mod(f(m, n) - f(m, max(1, n - 2)), 256) > Threshold
            h(f(m, n) + 1) = h(f(m, n) + 1) + 1;
            count = count + 1;
        end
    end
end

% Normalize gk to get k*
k_star = k / 256;

% Set a minimum value for u
umin = 0.001; % You can adjust this value

% Normalize count to get u
u = min(count / 256, umin);

% step7:Loop through each bin in the histogram
for bin = 1:256
    % Check condition and update histogram
    if b < bin && bin < w
        h(bin) = (1 - k_star) * u + k_star * h(bin);
    else
        alpha = 1; % You can adjust this value
        h(bin) = (1 / (1 + alpha)) * ((1 - k_star) * u + k_star * h(bin));
    end
end

% step8:obtain image from modified histogram

out=hist2image(h,M1,N1,b,w,alpha,I);
figure,imshow(out),title('Proposed Method Output')
figure,imhist(out(:,:,1));title(['Histogram of Proposed Method']);

% step9 : metrics calculation

% Convert images to grayscale if they are RGB
I = rgb2gray(I);
out = rgb2gray(out);
outHE=rgb2gray(outHE);

% Calculate AMBE (Absolute Mean Brightness Error)
ambe_HE = mean(mean(abs(double(I) - double(outHE))))
% Calculate Entropy (H) of low contrast image
entropyLowContrast = entropy(I)
% Calculate Entropy (H) of enhanced image
entropyEnhanced_HE = entropy(outHE)

% parameters for proposed method

% Calculate AMBE (Absolute Mean Brightness Error)
ambe = mean(mean(abs(double(I) - double(out))))
% Calculate Entropy (H) of low contrast image
entropyLowContrast = entropy(I)
% Calculate Entropy (H) of enhanced image
entropyEnhanced = entropy(out)





