function res_mser=mser_det(I,fg_mask,delta,allowed_area_change,min_object_size)

% VLFeat package needed
% http://www.vlfeat.org/download.html


%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       delta: MSER treshold step
%       allowed_area_change: MSER parameter
%       min_object_size:

%Output:
%       res_mser: image with binary points = detected cell centers

%Reference:
% [1] Matas, J., Chum, O., Urban, M., Pajdla, T.: Robust wide-baseline stereo from maximally stable extremal
% regions. Image and Vision Computing vol. 22(issue 10), 761–767 (2004). doi:10.1016/j.imavis.2004.02.006
% [2] Arteta, C., Lempitsky, V., Noble, J.A., Zisserman, A.: Learning to Detect Cells Using Non-overlapping Extremal
% Regions. Miccai (Figure 1), 1–8 (2012). doi:10.1007/978-3-642-33415-3



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%cell detection based on MSER - extract smalest region from overlaped ones










I=uint8(I*255);
% 
% regions = detectMSERFeatures(I,'ThresholdDelta',delta,'MaxAreaVariation',allowed_area_change,'RegionAreaRange',[min_object_size 999999]);
% 
% for k=1:length(regions)
%     pixely=regions(k).PixelList;
%     oblast=zeros(size(I));
%     for kk=1:size(pixely,1)
%         oblast(pixely(kk,2),pixely(kk,1))=1;
%     end
%     vsechny(:,:,k)=oblast;
% end
% 
% 



 [r,f]=vl_mser(I,'MinDiversity',0.2,...
        'MaxVariation',allowed_area_change,...
        'Delta',delta,...
        'MinArea', min_object_size/ numel(I),...
        'MaxArea',1);
    
    M1 = zeros(size(I)) ;
    for x=1:length(r)
        s = vl_erfill(I,r(x)) ;
        M1(s) = M1(s) + 1;
    end




suma=M1;


% suma=sum(vsechny,3);

m=imregionalmax(suma);

regions=m;
points=zeros(size(regions));
s = regionprops(regions>0,'centroid');
centroids = round(cat(1, s.Centroid));
for kp=1:size(centroids,1)
    points(centroids(kp,2),centroids(kp,1))=1;
end




res_mser=points.*fg_mask;
