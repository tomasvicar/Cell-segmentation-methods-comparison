function [detected]=hlog_zhang(a,gamma,fg_mask)


%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       gamma: normalization parameter - preferences of smaler (<2) or larger (>2) cells;


%Output:
%       detected: image with binary points = detected cell centers


%Reference:
% [1] Zhang, M., Wu, T., Bennett, K.M.: Small blob identification in medical images using regional features from
% optimum scale. IEEE Transactions on Biomedical Engineering vol. 62(issue 4), 1051–1062 (2015).
% doi:10.1109/TBME.2014.2360154

%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%Laplacian of Gaussian for cell detection - with optimal scale selection 
% based on Hesian computation 






sigma=1:1:35;


for k=1:length(sigma)
    
    
    sig=sigma(k);
    velikost= 2*ceil(3*sig)+1;
    h=(sig.^gamma)*fspecial('log', velikost, sig);
    
    
    pom=conv2_spec_symetric(a,h);
    log_all(:,:,k)=pom;
    
    I(:,:,k)=positive_def(pom);%.*fg_mask;
    pom=pom(I(:,:,k)==1);
    kvality(k)=mean(abs(pom(:)));
    
end

[~,ind]=max(kvality);

regions=I(:,:,ind);


points=zeros(size(regions));
s = regionprops(regions>0,'centroid');
centroids = round(cat(1, s.Centroid));
for kp=1:size(centroids,1)
    points(centroids(kp,2),centroids(kp,1))=1;
end
detected=points.*fg_mask;


end










function I=positive_def(a)
[dx,dy]=gradient(a);
[dxx,dxy]=gradient(dx);
[~,dyy]=gradient(dy);

I=(dxx>0).*((dxx.*dyy-dxy.*dxy)>0);


end





