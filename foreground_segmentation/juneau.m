function [ bw ] = juneau( I,window,t, min_object_size)

%Input:
%       I: input cell image  
%       window: range filter window size
%       t: treshold
%       min_object_size: smaler binary objects wil be removed
%Output:
%       bw: segmented foreground 

%Reference:
% [1] Juneau, P.-M., Garnier, A., Duchesne, C.: Selection and tuning of a fast and simple phase-contrast microscopy
% image segmentation algorithm for measuring myoblast growth kinetics in an automated manner. Microscopy
% and microanalysis 19(4), 855–866 (2013)



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz



NHOOD= true(window);

J = rangefilt(I, NHOOD);

bw =J>t;

bw =imfill(bw ,'holes');
bw =bwareafilt(bw ,[min_object_size 9999999999999]);