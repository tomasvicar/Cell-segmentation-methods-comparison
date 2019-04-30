function res_glog_kong=glog_kong(I,fg_mask,sigma0,alpha)

%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       sigma0: roughly minimal cell perimeters
%       alpha: alpha parameter

%Output:
%       res_glog_kong: image with binary points = detected cell centers


%Reference:
% [1] Kong, H., Akakin, H.C., Sarma, S.E.: A generalized laplacian of gaussian filter for blob detection and its
% applications. IEEE Transactions on Cybernetics 43(6), 1719–1733 (2013). doi:10.1109/TSMCB.2012.2228639



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%Laplacian of Gaussian for cell detection - generalized for
%eliptic shape cells



t=log(sigma0):0.2:3.5;
theta=pi/4;

sigma=exp(t);


for k=1:length(sigma)
    
    gamma=2;
    sig=sigma(k);
    filter_size= 2*ceil(3*sig)+1;
    h=(sig.^gamma)*fspecial('log', filter_size, sig);
    
    
    pom=conv2_spec_symetric(I,h);
    log_all(:,:,k)=-1.*pom;
    
    
end


l_g=max(log_all(:));
for k=1:length(sigma)
    s_i_l=min(min(log_all(:,:,k)));
    
    SS_i=(log_all(:,:,k)-s_i_l)/(l_g-s_i_l);
    zeta(k)=sum(sum(SS_i>0.6));
end
[~,best]=max(zeta);



sigma_x_min=sigma(max([best-3,1]));
sigma_x_max=sigma(min([best+3,length(sigma)]));

sigma_x=ceil(sigma_x_min):floor(sigma_x_max);

R=zeros(size(I));



Theta=linspace(0, pi - theta, round(pi / theta));
for sx=sigma_x
    
    sigma_y=sigma0:sx-1;
    
    for sy=sigma_y
        
        for theta=Theta

            
            nor=(1+log(sx)^alpha)*(1+log(sy)^alpha);
            
            [LoG]= glogkernel(sx, sy, theta);
            
            LoGn=nor*LoG;
            
            tmp=conv2_spec_symetric(I,LoGn);
            R=R-1.*tmp;
            
        end
        
        normalization=(1+log(sx)^alpha)*(1+log(sx)^alpha);
        
        [LoG]= glogkernel(sx, sx, theta);
        
        LoGn=normalization*LoG;
        
        tmp=conv2_spec_symetric(I,LoGn);
        R=R-1.*tmp;
        
    end
    
end


points=imregionalmax(R);
points=points.*fg_mask;
res_glog_kong=points;




end







function [LoG]= glogkernel(sigma_x, sigma_y, theta)

N =  ceil(2 * 3 * sigma_x);



[X, Y] =  meshgrid( linspace(0, N, N + 1) - N/2, linspace(0, N, N + 1) - N / 2);
a =  cos(theta) ^ 2 / (2 * sigma_x ^ 2) +  sin(theta) ^ 2 / (2 * sigma_y ^ 2);
b = - sin(2 * theta) / (4 * sigma_x ^ 2) + sin(2 * theta) / (4 * sigma_y ^ 2);
c =  sin(theta) ^ 2 / (2 * sigma_x ^ 2) + cos(theta) ^ 2 / (2 * sigma_y ^ 2);
D2Gxx = ((2*a*X + 2*b*Y).^2 - 2*a) .*  exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
D2Gyy = ((2*b*X + 2*c*Y).^2 - 2*c) .*  exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
Gaussian =  exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
LoG = (D2Gxx + D2Gyy) ./  sum(Gaussian(:));



end







