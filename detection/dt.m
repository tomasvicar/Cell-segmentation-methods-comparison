function res_dt=dt(I,fg_mask,min_object_size,min_hole_size,h,t)


%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       min_object_size: smaler binary objects will be removed
%       min_hole_size: smaler holes will be removed
%       h: h parameter of the h-minima transform
%       t: treshold value

%Output:
%       res_dt: image with binary points = detected cell centers

%Reference:
% [1] Thirusittampalam, K., Hossain, M.J., Ghita, O., Whelan, P.F.: A novel framework for cellular tracking and
% mitosis detection in dense phase contrast microscopy images. IEEE Journal of Biomedical and Health
% Informatics 17(3), 642–653 (2013). doi:10.1109/TITB.2012.2228663



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz

%cell detection based on h-minama and distance trasform of treshlded image







m=I>t;


m=~bwareafilt(~m,[round(min_hole_size)  inf]);
m = bwareafilt(m,[round(min_object_size),Inf]);

m=bwdist(m==0);
m=imhmax(m,h);

m=imregionalmax(m);



objects=m;
points=zeros(size(objects));
s = regionprops(objects>0,'centroid');
centroids = round(cat(1, s.Centroid));
for kp=1:size(centroids,1)
    points(centroids(kp,2),centroids(kp,1))=1;
end

res_dt=points.*fg_mask;

