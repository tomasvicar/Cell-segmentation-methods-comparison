% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.



% S = auto_edge_detection_ksdensity(I, min_cell_size, lower_hole_size_bound, upper_hole_size_bound, manual_finetune)
% 
% Inputs
%     -I: image to be segmented using automatic gradient based edge detection.
%     -min_cell_size: smallest object in segmented image that will be kept. Smaller objects will be removed.
%     -upper_hole_size_bound: holes larger than this in the segmented image will be filled in.
%     -manual_finetune: 
% 
% Outputs
%     -S: segmented image (logical)

function S = EGT_Segmentation(I, min_cell_size, upper_hole_size_bound, manual_finetune)
if ~exist('min_cell_size','var')
    min_cell_size = 1;
end
if ~exist('upper_hole_size_bound','var')
    upper_hole_size_bound = inf;
end
if ~exist('manual_finetune','var')
    manual_finetune = 0;
end

% this controls how far each increment of manual_finetune moves the percentile threshold
greedy_step = 1; 

% Compute image gradient and percentiles
S = imgradient(single(I));
S1 = nonzeros(S);

% Compute gradient image histogram
ratio = (max(S1)-min(S1))/1000;
% factor = ceil(numel(S1)/25000000); % set the factor so limit the number of pixel to 25 million
% [hist_data,~] = hist(S1(1:factor:end),min(S1):ratio:max(S1));
[hist_data,~] = hist(S1(:),min(S1):ratio:max(S1));


% get the mode and the corresponding frequency value
[~,hist_mode_loc] = sort(hist_data, 'descend');
hist_mode_loc = round(mean(hist_mode_loc(1:3))); % take the average of the first 3 peaks
% normalize the histogram counts between 0-1
temp_hist = hist_data/sum(hist_data)*100;
% compute lower bound
lower_bound = 3*hist_mode_loc;
if lower_bound > numel(temp_hist)
    warning('lower bound set to end of list.');
    lower_bound = numel(temp_hist);
end

% ensure that 75% of the pixels have been taken
% c = cumsum(temp_hist);
% idx = find(c>95,1);
norm_hist = temp_hist/max(temp_hist);
idx = find(norm_hist(hist_mode_loc:end)<0.05,1) + hist_mode_loc - 1;

upper_bound = max(idx, 18*hist_mode_loc);

% Compute the density metric
if upper_bound > numel(temp_hist)
    warning('upper bound set to end of list.');
    upper_bound = numel(temp_hist); 
end
density_metric = sum(temp_hist(lower_bound:upper_bound));

% Fit a line between the 80th and the 40th percentiles from the plot above
saturation1 = 3;
saturation2 = 42;
a = (95 - 40) / (saturation1 - saturation2);
b = 95 - a*saturation1;

% Compute gradient threshold
prct_value = round(a*density_metric + b);
if prct_value > 98;
    prct_value = 98;
end
if prct_value < 25;
    prct_value = 25;
end
% decrease or increase by a multiple of 5 percentile the manual input
prct_value = prct_value - greedy_step*manual_finetune; 
if prct_value > 100, prct_value = 100; end
if prct_value < 1, prct_value = 1; end
	
% threshold = prctile(S1,prct_value);
prct_value = prct_value/100;
threshold = percentile_computation(S1,prct_value);


% Threshold the gradient image and perform some cleaning with morphological operations
S = S > threshold;
S = fill_holes(S, upper_hole_size_bound);
S = imerode(S,strel('disk',1));  % removes the border pixels from the gradient
S = bwareaopen(S,min_cell_size,8);

end





function T = percentile_computation(A,p)
% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

assert(all(p >= 0 & p <= 1), 'Percentiles must be between zero and one inclusive.');

B = A(~isnan(A)); % ignore NaNs
if isempty(B)
    T = NaN(size(p));
    return;
end

B = sort(B(:));

indx = round(p*length(B)+1);

indx(indx<1) = 1;
indx(indx>length(B)) = length(B);

T = B(indx);

T = reshape(T, size(p));

end







function S = fill_holes(S, upper_bound)

% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.


S = logical(S);

if isinf(upper_bound)
    % the user has selected no upper limit (inf) for the fill holes size
    CC = bwconncomp(~S, 4);
    szs = zeros(numel(CC.PixelIdxList),1);
    for i = 1:numel(CC.PixelIdxList)
        szs(i,1) = numel(CC.PixelIdxList{i});
    end
    upper_bound = max(szs(:)) - 1;
    if isempty(upper_bound)
        upper_bound = 0;
    end
end

% BWu holds the mask of pixels where the holes larger than upper_bound have been filled
BWu = bwareaopen(~S, upper_bound, 4);
BWu = (~S) - BWu;

% remove the hole pixels from S that are in either the upper or lower mask
S(BWu>0) = 1;

end








