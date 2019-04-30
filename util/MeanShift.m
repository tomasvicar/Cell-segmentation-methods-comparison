function [clustCent] = MeanShift(dataPts,bandWidth)

%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz




stop_thresh = 1e-3;
num_points=size(dataPts,2);
num_clust=0;

permuted_ind=randperm(num_points);
citac=0;
for ind=permuted_ind
    new=dataPts(:,ind);
    citac=citac+1;
%     disp([num2str(citac) '/' num2str(num_points)])
    moving=1;
    while moving
        old=new;
        sq_dist =  sum(bsxfun(@minus,new,dataPts).^2);
        
        used_points = find(sq_dist < bandWidth^2);
        
        new=mean(dataPts(:,used_points),2);
        
        
        
        
        if abs(new-old) < stop_thresh
            moving=0;
            mergeWith=0;
            for c = 1:num_clust
                if norm(new-clustCent(:,c)) < bandWidth/2 
                    
                    mergeWith = c;
                    
                    break;
                end
                
            end
            if mergeWith > 0 
                clustCent(:,mergeWith) = new;
            else 
                num_clust = num_clust+1; 
                clustCent(:,num_clust) = new; 
            end
            
         
        end
        
        
        
        
        
        
    end
    
    
    
end





%
%
% numClust = 0;
% stopThresh = 1e-3*bandWidth; % when mean has converged
% clustCent = []; % center of clust
%
% numPts=size(dataPts,2); % number of points to posibaly use as initilization points
% clusterVotes = zeros(1,numPts,'uint16'); % used to resolve conflicts on cluster membership
%
%
% kmean = @(x,dis) mean(x,2);
%
%
% usedPoints=false(1,numPts);
%
% numOfAllPts=numPts;
%
%
% for k=1:numOfAllPts
%
%     initPtInds=find(~usedPoints);
%     numPts=lenngth(initPtInds);
%     tempInd = ceil( (numPts-1e-6)*rand);
%     stInd = initPtInds(tempInd); % use this point as start of mean
%     myMean = dataPts(:,stInd);  % intilize mean to this points location
%     myMembers = []; % points that will get added to this cluster
%     thisClusterVotes = zeros(1,numPts,'uint16'); % used to resolve conflicts on cluster membership
%
%     while true %loop untill convergence
%         sqDistToAll = sum(bsxfun(@minus,myMean,dataPts).^2); % dist squared from mean to all points still active
%
%         inInds = find(sqDistToAll < bandWidth^2); % points within bandWidth
%         thisClusterVotes(inInds) = thisClusterVotes(inInds)+1; % add a vote for all the in points belonging to this cluster
%
%         myOldMean = myMean; % save the old mean
%         myMean = kmean(dataPts(:,inInds),sqrt(sqDistToAll(inInds))); % compute the new mean
%         myMembers = [myMembers inInds]; % add any point within bandWidth to the cluster
%         usedPoints(myMembers) = true; % mark that these points have been visited
%
%
%
%         %**** if mean doesn't move much stop this cluster ***
%         if norm(myMean-myOldMean) < stopThresh
%             %check for merge posibilities
%             mergeWith = 0;
%             for cN = 1:numClust
%                 distToOther = norm(myMean-clustCent(:,cN)); % distance to old clust max
%                 if distToOther < bandWidth/2 % if its within bandwidth/2 merge new and old
%                     mergeWith = cN;
%                     break;
%                 end
%             end
%
%             if mergeWith > 0 % something to merge
%                 nc = numel(myMembers); % current cluster's member number
%                 no = numel(clustMembsCell{mergeWith}); % old cluster's member number
%                 nw = [nc;no]/(nc+no); % weights for merging mean
%                 clustMembsCell{mergeWith} = unique([clustMembsCell{mergeWith},myMembers]);   %record which points inside
%                 clustCent(:,mergeWith) = myMean*nw(1) + myOldMean*nw(2);
%                 clusterVotes(mergeWith,:) = clusterVotes(mergeWith,:) + thisClusterVotes;    %add these votes to the merged cluster
%             else % it's a new cluster
%                 numClust = numClust+1; %increment clusters
%                 clustCent(:,numClust) = myMean; %record the mean
%                 clustMembsCell{numClust} = myMembers; %store my members
%                 clusterVotes(numClust,:) = thisClusterVotes; % creates a new vote
%             end
%
%             break;
%         end
%
%     end
%
% end
%
%
%
% end
