function S=gfrst(a,fg_mask,rrange,t,kr,alpha)

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
% [1] Bahlmann, C.: Fast radial symmetry detection under affine transformations. In: Proceedings of the 2012 IEEE
% Conference on Computer Vision and Pattern Recognition (CVPR). CVPR ’12, pp. 932–939. IEEE Computer
% Society, Washington, DC, USA (2012). http://dl.acm.org/citation.cfm?id=2354409.2354741



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% % vicar@feec.vutbr.cz


%fast radial symmetry tranform for cell detection - generalized for
%eliptic shape cells



aa = padarray(a,[max(rrange)*2+10 max(rrange)*2+10]);


S=zeros([size(aa) length(kr) length(alpha)]);




[gxx,gyy] = gradient(a);
ag=sqrt(gxx.^2+gyy.^2);

ag = padarray(ag,[max(rrange)*2+10 max(rrange)*2+10]);




above_tresh=ag>0;
[x,y]=find(above_tresh);
ag=ag(above_tresh);



theta_all=(0:5)*pi/6;
a_all=[6, 8, 10, 12, 14, 16];
b_all =[4, 6, 8];

for a=a_all
    for b=b_all
        
        for theta=theta_all
            
            
            
            G=[cos(theta),-sin(theta);sin(theta),cos(theta)]*[a 0; 0 b];
            M=[0 1;-1 0];
            
            
            
            g=(G*M*G^(-1)*M^(-1))*[gxx(:)';gyy(:)'];
            
            gy=reshape(g(1,:),size(gxx));
            gx=reshape(g(2,:),size(gyy));
            
            agg=sqrt(gx.^2+gy.^2);
            
            ngx=gx./(agg+eps);
            
            ngy=gy./(agg+eps);
            
            
            ngx = padarray(ngx,[max(rrange)*2+10 max(rrange)*2+10]);
            
            ngy = padarray(ngy,[max(rrange)*2+10 max(rrange)*2+10]);
            
            
            
            ngx=ngx(above_tresh);
            ngy=ngy(above_tresh);
            
            
            
            
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
                
                %     toc
                
            end
        end
    end
end
S=S(max(rrange)*2+10+1:end-max(rrange)*2-10,max(rrange)*2+10+1:end-max(rrange)*2-10,:,:);




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


