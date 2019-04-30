function y=conv2_spec_symetric(x,h)

%Input:
%       x: input image
%       h: convolution mask


%Output:
%       y: result image


%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@fee.vutbr.cz


%fast spectral domain convolution




n=floor(size(h,1)/2);
m=floor(size(h,2)/2);
x = padarray(x,[n m],'symmetric');

y=real(ifft2(fft2(x,size(x,1),size(x,2)).*fft2(h,size(x,1),size(x,2))));
y=circshift(y,floor(-1*[size(h)]/2));

y=y(n+1:end-n,m+1:end-m);


end
