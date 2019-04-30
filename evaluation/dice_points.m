function dice=dice_points(gt_mask,detection)


%Input:
%       gt_mask: ground true segmentation mask
%       detection: cell detection results - bw image, where cells are binary points


%Output:
%       dice: dice coeficient


%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz




TP=0;
FN=0;
FP=0;

gt_mask=gt_mask==1;
detection=detection==1;

for kk=1:size(gt_mask,3)
    mm=bwareafilt(gt_mask(:,:,kk),[50,999999]);
    l=bwlabel(mm,4);
    det=detection(:,:,kk);

for k=1:max(l(:))
    cell=k==l;
    if sum(sum(cell.*det))>0
        TP=TP+1;
        FP=FP+sum(sum(cell.*det))-1;
    else
        FN=FN+1;
    end
end
FP=FP+sum(sum(det.*(mm==0)));
end

dice=2*TP/(2*TP+FP+FN);
precision=TP/(TP+FP);
recall=TP/(TP+FN);
jaccard=dice/(2-dice);
% accuracy=(TP+TN)/(TP+TN+FP+FN);

end