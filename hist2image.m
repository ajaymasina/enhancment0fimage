
function out_RGB=hist2image(temp,M,N,b,w,alpha,I)
mu=alpha;

[Y,U,V] = rgb2yuv(I(:,:,1),I(:,:,2),I(:,:,3));
Yn = double(Y);
temp=0;

% extract an input histogram vector h
h = zeros(256,1);
for j=1:M
    for i=1:N
        temp = Yn(j,i) + 1;
        h(temp,1) = h(temp,1) + 1;
    end
end
clear temp
mu=[2];
m = LHM(h, mu);
% Convex optimization
D = inv(tril(ones(256,256)));
[y, ~, ~] = histmodification(m, h, 1); x = D\y;
% write output image
out_Y= zeros(M,N);
for j=1:M
    for i=1:N
        out_Y(j,i) = round(x(Yn(j,i)+1));
    end
end

out_RGB = yuv2rgb(uint8(out_Y),U,V);
mm=0.01*(b/2*w);
out_RGB=out_RGB*mm;
