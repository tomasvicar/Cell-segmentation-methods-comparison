function res_glog_xu=glog_xu(I,fg_mask,sigma,mean_shift_window)


%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       sigma: cell perimeters range rmin/sqrt(2):rmax/sqrt(2);
%       mean_shift_window:

%Output:
%       res_glog_xu: image with binary points = detected cell centers


%Reference:
% [1] Xu, H., Lu, C., Berendt, R., Jha, N., Mandal, M.: Automatic nuclei detection based on generalized laplacian of
% gaussian filters. IEEE journal of biomedical and health informatics 21(3), 826–837 (2017)


%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%Laplacian of Gaussian for cell detection - generalized for
%eliptic shape cells



theta=pi/9;
sigma_x=sigma;
a=I;
Theta=linspace(0, pi - theta, round(pi / theta));
counter=0;
total_log=zeros(201,201,length(Theta)+1);
for theta=Theta
    counter=counter+1;
    
    for sx=sigma_x
        %
        sigma_y=sigma(1):sx;
        %
        for sy=sigma_y
            
            nor=sx*sy;
            
            [LoG]= glogkernel(sx, sy, theta);
            
            LoGn=nor*LoG;
            
            total_log(:,:,counter)=total_log(:,:,counter)+LoGn;
            
            
        end
    end
end
counter=counter+1;

for sx=sigma_x

    nor=sx*sy;
    [LoG]= glogkernel(sx, sx, theta);
    
    LoGn=nor*LoG;
    
    total_log(:,:,counter)=total_log(:,:,counter)+LoGn;
end



detected=zeros(size(a));
for k=1:size(total_log,3)
    pom=conv2_spec_symetric(a, total_log(:,:,k));
    detected=detected+imregionalmin(pom);
end
detected=detected.*fg_mask;


[x,y]=find(detected);
data=[x,y]';
[clustCent] = MeanShift(data,mean_shift_window);
clustCent=round(clustCent)';


points=zeros(size(a));
for k=1:size(clustCent,1)
    points(clustCent(k,1),clustCent(k,2))=1;
end

res_glog_xu=points;


end









function [LoG]= glogkernel(sigma_x, sigma_y, theta)


N=200;

[X, Y] =  meshgrid( linspace(0, N, N + 1) - N/2, linspace(0, N, N + 1) - N / 2);
a =  cos(theta) ^ 2 / (2 * sigma_x ^ 2) +  sin(theta) ^ 2 / (2 * sigma_y ^ 2);
b = - sin(2 * theta) / (4 * sigma_x ^ 2) + sin(2 * theta) / (4 * sigma_y ^ 2);
c =  sin(theta) ^ 2 / (2 * sigma_x ^ 2) + cos(theta) ^ 2 / (2 * sigma_y ^ 2);
D2Gxx = ((2*a*X + 2*b*Y).^2 - 2*a) .*  exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
D2Gyy = ((2*b*X + 2*c*Y).^2 - 2*c) .*  exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
Gaussian =  exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
LoG = (D2Gxx + D2Gyy) ./  sum(Gaussian(:));



end











