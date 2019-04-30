function t=poisson_tresh(I)


%Input:
%       I: input image
%Output:
%       bw: segmented foreground 

%Reference:
% [1] Al-Kofahi, Y., Lassoued, W., Lee, W., Roysam, B.: Improved automatic detection and segmentation of cell
% nuclei in histopathology images. IEEE Transactions on Biomedical Engineering 57(4), 841–852 (2010).
% doi:10.1109/TBME.2009.2035102



%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz

%minimum error tresholding for poissson distribution



n_bins=128;


I_norm=mat2gray(I);


edges=linspace(0,1,n_bins);
N= histcounts(I_norm(:),edges, 'Normalization', 'probability');


for t=1:length(N)
    P0(t)=sum(N(1:t));
    P1(t)=sum(N(t+1:end));
    
    mu0(t)=sum(((1:t)-1).*N(1:t))/P0(t);
    mu1(t)=sum(((t+1:numel(N))-1).*N(t+1:end))/P1(t);
    
end


mu=sum(((1: numel(N))-1).*N);

tt=mu-(P0.*(log(P0)+mu0.*log(mu0)))-(P1.*(log(P1)+mu1.*log(mu1)));


[~,t]=min(tt);

t=edges(t);
t=t*(max(I(:))-min(I(:)))+min(I(:));



drawnow