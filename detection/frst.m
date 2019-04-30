function [S] = frst(I,fg_mask,rrange,t,kr,alpha)

%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       rrrange: range of cell perimeters
%       t: treshold of votting image
%       kr: k parameter
%       alpha: alpha parameter

%Output:
%       S: image with binary points = detected cell centers


%Reference:
% [1] Loy, G., Zelinsky, A.: Fast radial symmetry for detecting points of interest. IEEE Transactions on Pattern
% Analysis and Machine Intelligence 25(8), 959–973 (2003). doi:10.1109/TPAMI.2003.1217601



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% % vicar@feec.vutbr.cz


%fast radial symmetry tranform for cell detection







[gy,gx] = gradient(I);

ag=sqrt(gx.^2+gy.^2);

ngx=gx./(ag+eps);

ngy=gy./(ag+eps);


aa = padarray(I,[max(rrange)+5 max(rrange)+5]);

ag = padarray(ag,[max(rrange)+5 max(rrange)+5]);

ngx = padarray(ngx,[max(rrange)+5 max(rrange)+5]);

ngy = padarray(ngy,[max(rrange)+5 max(rrange)+5]);

agg=ag;


above_tresh=ag>0;

[x,y]=find(above_tresh);

ngx=ngx(above_tresh);
ngy=ngy(above_tresh);

ag=ag(above_tresh);

S=zeros([size(aa) length(kr) length(alpha)]);


k=0;
for r = rrange
    k=k+1;
    
    Mr=zeros([size(aa)]);
    Or=zeros([size(aa)]);
    
    
    
    
    px=x+round(ngx*r);
    py=y+round(ngy*r);
    
    
    pom=sub2ind(size(Or), px,py);
    
    
    Or =accumarray( pom, ones(1,length(pom)),[numel(Or) 1]);
    Or=reshape(Or,size(aa));
    
    
    Mr =accumarray( pom,ag,[numel(Mr) 1]);
    Mr=reshape(Mr,size(aa));
    
    
    
    
    Or_hat=repmat(Or,[1 1 length(kr)]);
    clear krr
    krr(1,1,1:length(kr))=kr;
    krr=repmat(krr,[size(Or) 1]);
    Or_hat(Or_hat>krr)=krr(Or_hat>krr);
    Mrr=repmat(Mr,[1 1 length(kr)]);
    
    for qq=1:length(alpha)
        Fr=(Mrr./krr).*(Or_hat./krr).^(alpha(qq));
        S(:,:,:,qq)=S(:,:,:,qq)+imgaussfilt(Fr,0.5*r);
    end
    
    
    
end
S=S(max(rrange)+5+1:end-max(rrange)-5,max(rrange)+1+5:end-max(rrange)-5,:,:);







objects=imregionalmax(S);

objects(S<t)=0;
points=zeros(size(objects));
s = regionprops(objects>0,'centroid');
centroids = round(cat(1, s.Centroid));
for kp=1:size(centroids,1)
    points(centroids(kp,2),centroids(kp,1))=1;
end
points(centroids(kp,2),centroids(kp,1))=1;


points=points.*fg_mask;




S=points;




end





