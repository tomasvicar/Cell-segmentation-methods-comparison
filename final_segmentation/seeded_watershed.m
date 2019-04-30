function final_seg=seeded_watershed(I,points,fg_mask)

%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       points: image with binary points = detected cell centers


%Output:
%       final_seg: binary image with final single cell segmetnation

%Reference:
% [1] Parvati, K., Rao, P., Mariya Das, M.: Image segmentation using gray-scale morphology and marker-controlled
% watershed transformation. Discrete Dynamics in Nature and Society 2008 (2008)



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz

%marker-controled watershed segmentation







imp = imimposemin(-I, points);


final_seg=double(watershed(imp)>0).*double(fg_mask);

final_seg = -imimposemin(-final_seg, points)>0;

