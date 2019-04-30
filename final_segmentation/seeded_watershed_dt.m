function final_seg=seeded_watershed_dt(I,points,fg_mask)


%Input:
%       I: input image
%       fg_mask: binary foreground mask
%       points: image with binary points = detected cell centers


%Output:
%       final_seg: binary image with final single cell segmetnation

%Reference:
% [1] Parvati, K., Rao, P., Mariya Das, M.: Image segmentation using gray-scale morphology and marker-controlled
% watershed transformation. Discrete Dynamics in Nature and Society 2008 (2008)
% [2] Ikonen, L., Toivanen, P.: Shortest routes on varying height surfaces using gray-level distance transforms. Image
% and Vision Computing 23(2), 133–141 (2005)


%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%marker-controled watershed segmentation of distance transform image (geodesic distance)




dt=bwdistgeodesic(fg_mask>0,points>0);
dt(isnan(dt))=99999;
fg_mask(dt==Inf)=0;

final_seg=double(watershed(dt)>0).*double(fg_mask);




