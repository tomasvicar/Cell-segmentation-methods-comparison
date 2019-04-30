function points=mlog_peng(I,fg_mask,sigma,gamma,t)

%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       sigma: cell perimeters range rmin/sqrt(2):rmax/sqrt(2);
%       gamma: normalization parameter - preferences of smaler (<2) or larger (>2) cells;
%       t: treshold

%Output:
%       points: image with binary points = detected cell centers


%Reference:
% [1] Peng, H., Zhou, X., Li, F., Xia, X., Wong, S.T.C.: Integrating multi-scale blob/curvilinear detector techniques
% and multi-level sets for automated segmentation of stem cell images. In: 2009 IEEE International Symposium
% on Biomedical Imaging: From Nano to Macro, pp. 1362–1365. IEEE, ??? (2009)

%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%Laplacian of Gaussian for cell detection 



log_map=zeros([size(I) length(sigma)]);
for k=1:length(sigma)
    
    
    sig=sigma(k);
    filter_size= 2*ceil(3*sig)+1;
    h=(sig.^gamma)*fspecial('log', filter_size, sig);
    
    
    pom=conv2_spec_symetric(I,h);
    log_map(:,:,k)=pom;
    
    
end
log_map=min(log_map,[],3);


bw = log_map<t;



points=zeros(size(bw));
s = regionprops(bw>0,'centroid');
centroids = round(cat(1, s.Centroid));
for kp=1:size(centroids,1)
    points(centroids(kp,2),centroids(kp,1))=1;
end

points=points.*fg_mask;





end














