function dice=dice_foreground(gt_mask,foreground)



%Input:
%       gt_mask: ground true mask
%       foreground: binary foreground segmentation mask


%Output:
%       dice: dice coeficient


%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz














gt_mask=gt_mask>0;
foreground=foreground>0;
TP=sum(sum(sum(gt_mask.*foreground)));

TN=sum(sum(sum((gt_mask==0).*(foreground==0))));

FN=sum(sum(sum((gt_mask==1).*(foreground==0))));

FP=sum(sum(sum((gt_mask==0).*(foreground==1))));



dice=2*TP/(2*TP+FP+FN);
precision=TP/(TP+FP);
recall=TP/(TP+FN);
jaccard=dice/(2-dice);
accuracy=(TP+TN)/(TP+TN+FP+FN);



end