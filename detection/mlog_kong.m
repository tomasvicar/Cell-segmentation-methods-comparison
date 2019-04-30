function res_mlog_kong=mlog_kong(I,fg_mask,sigma,gamma,t,max_overlap)


%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       sigma: cell perimeters range rmin/sqrt(2):rmax/sqrt(2);
%       gamma: normalization parameter - preferences of smaler (<2) or larger (>2) cells;
%       t: treshold
%       max_overlap: alowd ovetlap for regions pruning

%Output:
%       res_mlog_kong: image with binary points = detected cell centers


%Reference:
% [1] Kong, H., Akakin, H.C., Sarma, S.E.: A generalized laplacian of gaussian filter for blob detection and its
% applications. IEEE Transactions on Cybernetics 43(6), 1719–1733 (2013). doi:10.1109/TSMCB.2012.2228639

%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%Laplacian of Gaussian for cell detection - method without generalization
%from [1]



for k=1:length(sigma)
    
    
    sig=sigma(k);
    filter_size= 2*ceil(3*sig)+1;
    h=(sig.^gamma)*fspecial('log', filter_size, sig);
    
    
    tmp=conv2_spec_symetric(I,h);
    log_all(:,:,k)=tmp;
    
    
    
    
end


points=imregionalmin(log_all);
points=points.*repmat(fg_mask,[1 1 size(points,3)]);



points(log_all>t)=0;





[r,c,v] = ind2sub(size(points),find(points));
for kp=1:length(r)
    for kpp=kp+1:length(r)
        area=itersect_circel(r(kp), c(kp), sqrt(2)*sigma(v(kp)), r(kpp), c(kpp), sqrt(2)*sigma(v(kpp)));
        area1=2*pi*sqrt(2)*sigma(v(kp));
        area2=2*pi*sqrt(2)*sigma(v(kpp));
        if area>(max_overlap*min([area1,area2]))
            if area1<area2
                points(r(kp), c(kp),v(kp))=0;
            else
                points(r(kpp), c(kpp),v(kpp))=0;
            end
        end
        
        
    end
end


points2=sum(points,3);



points=zeros(size(points2));
s = regionprops(points2>0,'centroid');
centroids = round(cat(1, s.Centroid));
for kp=1:size(centroids,1)
    points(centroids(kp,2),centroids(kp,1))=1;
end




res_mlog_kong=points;




end




















function area=itersect_circel(x0, y0, r0, x1, y1, r1)

rr0 = r0 * r0;
rr1 = r1 * r1;
d = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));


if (d > r1 + r0)
    
    area=0;
    
elseif (d <= abs(r0 - r1) && r0 >= r1)
    
    
    area= pi* rr1;
    
    
    
    
    
elseif (d <= abs(r0 - r1) && r0 < r1)
    
    
    area=pi * rr0;
    
    
    
    
else
    
    phi = (acos((rr0 + (d * d) - rr1) / (2 * r0 * d))) * 2;
    theta = (acos((rr1 + (d * d) - rr0) / (2 * r1 * d))) * 2;
    area1 = 0.5 * theta * rr1 - 0.5 * rr1 * sin(theta);
    area2 = 0.5 * phi * rr0 - 0.5 * rr0 * sin(phi);
    
    
    
    area= area1 + area2;
end

end




