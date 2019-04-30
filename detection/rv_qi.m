function res_rv_qi=rv_qi(I,fg_mask,rrange,sigma1,w,sigma2)

%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       rrrange: range of cell perimeters
%       sigma1: Gaussian filter sigma
%       w: meanshift window size
%       sigma2: vote kernel sigma

%Output:
%       res_rv_qi: image with binary points = detected cell centers


%Reference:
% [1] Qi, X., Xing, F., Foran, D.J., Yang, L.: Robust segmentation of overlapping cells in histopathology specimens
% using parallel seed detection and repulsive level set. IEEE Transactions on Biomedical Engineering 59(3),
% 754–765 (2012)



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%single-pass radial votting cell detection



rmin=rrange(1);
rmax=rrange(end);

Delta=30;


S=0;

im=I;

rmean=round(mean([rmin rmax]));


im=imgaussfilt(im,sigma1);
[Gmag, Gdir] = imgradient(im,'intermediate');



imm=im;



im=padarray(im,[rmax+5, rmax+5]);
Gmag = padarray(Gmag,[rmax+5, rmax+5]);
Gdir = padarray(Gdir,[rmax+5, rmax+5]);


result=zeros(size(Gmag));

above_tresh=Gmag>S;

[x,y]=find(above_tresh);

Gmag=Gmag(above_tresh);

Gdir=Gdir(above_tresh);




a=im;
[grid_x,grid_y]=meshgrid(linspace(-size(a,2)/2,size(a,2)/2,size(a,2)),linspace(-size(a,1)/2,size(a,1)/2,size(a,1)));
cone=sqrt(grid_x.^2+grid_y.^2);
anglet=angle(grid_x - 1i*grid_y).*180/pi;
circle=(cone>rmin).*(cone<rmax);
kernel=(anglet>(-Delta)).*(anglet<(Delta)).*circle;
kernel=kernel>0;




% shift=round([-size(a,2)/2,-size(a,1)/2]);
result=zeros(size(a));
gg = fspecial('gaussian',size(a), sigma2);



Gdir=round(Gdir/3).*3;

 directions=-180:3:180;
 parfor pomk=1:length(directions)
    direction=directions(pomk) ;
    vykus=imrotate(kernel,direction,'crop');

    g=circshift(gg,round(rmean*[-sin(direction/180*pi),cos(direction/180*pi)]));
    vykus=vykus.*g;

    
    aktualni=Gdir==direction;
    px=x(aktualni);
    py=y(aktualni);
    Gmagg=Gmag(aktualni);

    vykuss=zeros(size(vykus));

    vykuss(sub2ind(size(vykus), px, py))=Gmagg;
    vykuss=conv2_spec_symetric(vykuss,vykus);
    result(:, :,pomk) =vykuss;
    
end
result=sum(result, 3);
result=result(rmax+5+1:end-rmax-5,rmax+1+5:end-rmax-5);
mmax=max(result(:));

% imshowpair(imm,vysledek)
x=[];
y=[];
for R=0.2:0.1:0.9
    tmp=(result>(R*mmax)).*fg_mask;
    [xx,yy]=find(tmp);
    x=[x;xx];
    y=[y;yy];
end




data=[x,y]';
[clustCent] = MeanShift(data,w);
clustCent=round(clustCent)';
points=zeros(size(imm));
for kk=1:size(clustCent,1)
points(clustCent(kk,1),clustCent(kk,2))=1;
end

res_rv_qi=points;

end




