function [dice]=dice_final_segmentation(gt_mask,segmentation)


%Input:
%       gt_mask: ground true mask
%       segmentation: binary single cell segmentation mask


%Output:
%       dice: mean dice coeficient


%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz


%Reference:
% [1] Ma¡ska, M., Ulman, V., Svoboda, D., Matula, P., Matula, P., Ederra, C., Urbiola, A., Espa?na, T., Venkatesan,
% S., Balak, D.M.W., Karas, P., Bolckov´a, T., Streitov´a, M., Carthel, C., Coraluppi, S., Harder, N., Rohr, K., ¡
% Magnusson, K.E.G., Jald´en, J., Blau, H.M., Dzyubachyk, O., K¡r´?¡zek, P., Hagen, G.M., Pastor-Escuredo, D.,
% Jimenez-Carretero, D., Ledesma-Carbayo, M.J., Mu?noz-Barrutia, A., Meijering, E., Kozubek, M.,
% Ortiz-de-Solorzano, C.: A benchmark for comparison of cell tracking algorithms. Bioinformatics vol. 30(issue
% 11), 1609–1617 (2014-6-1). doi:10.1093/bioinformatics/btu080


%mean dice coeficient - computed as in [1] with dice instead of jaccard














gt_mask=gt_mask>0;

counter=0;
for kk=1:size(gt_mask,3)
    
    YY=segmentation(:,:,kk);
    ll=bwlabel(YY);
    l=bwlabel(gt_mask(:,:,kk),4);
    for k=1:max(l(:))
        counter=counter+1;
        b=k==l;
        bb=ll(b);
        qq=unique(bb(find(bb)));
        dice(counter)=0;

        for q=qq'
            cell=(ll==q);
            if (sum(b(:))*0.5)<sum(sum((b&cell)))
                ll(cell)=0;
                dice(counter)=2*sum(sum((b&cell)))/(sum(cell(:))+sum(b(:)));
            end
                
        end

        
    end
    segmentation(:,:,kk)=YY;



end


dice=mean(dice);
