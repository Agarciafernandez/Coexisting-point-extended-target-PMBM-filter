function [tracks,table,wAssoc] = trackOrientedPMBPointExtended(tracks,table,wAssoc)

%number of tracks
n_tt = length(tracks);
if n_tt == 0 || length(wAssoc)==1
    return
else
    %Normalised wAssoc
    wAssocNorm=exp(wAssoc);
    wAssocNorm=log(wAssocNorm/sum(wAssocNorm));
    %merging local hypotheses in the same track
    for i = 1:n_tt
        %number of local hypotheses
        nb = length(tracks{i});
        w_merge = zeros(nb,1); %Angel: I use this one for GGIW merging (as was done in Yuxuan's code.
        GGIW_to_merge = [];
        r_merge=0;
        c_merge=0;
        w_merge_point=0;
        sum_weight_r=0;
        x_merged=zeros(4,1);
        P_merged=zeros(4);
        for j = 1:nb
            %add global hypotheses weights
                        
            [~,w_merge(j)] = normalizeLogWeights(wAssocNorm(table(:,i) == j));
            
            w_merge_j_exp=exp(w_merge(j));
            sum_weight_r=sum_weight_r+w_merge_j_exp*tracks{i}(j).Bern.r;
            
            if(isfield(tracks{i}(j).Bern,'Point') && tracks{i}(j).Bern.Point.c<1)
                
                w_merge_point_j=w_merge_j_exp*tracks{i}(j).Bern.r*tracks{i}(j).Bern.Point.c;
                xp=tracks{i}(j).Bern.Point.xp;
                Pp=tracks{i}(j).Bern.Point.Pp;
                x_merged=x_merged+w_merge_point_j*xp;
                P_merged=P_merged+w_merge_point_j*(xp*xp'+Pp);
                w_merge_point=w_merge_point+w_merge_point_j;
                %c_merge=c_merge+w_merge_j_exp*tracks{i}(j).Bern.Point.c;
                w_merge(j) = w_merge(j) + log(tracks{i}(j).Bern.r*(1-tracks{i}(j).Bern.Point.c)); %The merging of GGIW, must account for 1-c
            else
                %This hypothesis corresponds to an extended target
                w_merge(j) = w_merge(j) + log(tracks{i}(j).Bern.r);
            end
            
            
            r_merge=r_merge+w_merge_j_exp*tracks{i}(j).Bern.r;         
            %multiply with the existence probability (summation in logarithm)
            GGIW_to_merge = [GGIW_to_merge;tracks{i}(j).Bern.GGIW];
            
        end
        
        if(w_merge_point>0)
            %This means that at least one hypothesis can be a point
            %target so we have the Point field.
            x_merged=x_merged/w_merge_point;
            P_merged=P_merged/w_merge_point;
            P_merged=P_merged-(x_merged*x_merged');
            c_merge=w_merge_point/sum_weight_r;
            tracks{i}(1).Bern.Point.c = c_merge;
            tracks{i}(1).Bern.Point.xp=x_merged;
            tracks{i}(1).Bern.Point.Pp=P_merged;
        end
        
        if(~isreal(w_merge))
            display('error')
        end
        
        if(sum(w_merge==-Inf)<length(w_merge))
            %merge all GGIWs and existence probability
            [r_hat,GGIW_hat] = GGIW_merge_wrap(w_merge,GGIW_to_merge);
        else
            %The GGIW component has zero weight, but we must
            %include this field
            GGIW_hat=GGIW_to_merge(1);
        end
        
        tracks{i}(1).Bern.r=r_merge;
        tracks{i}(1).Bern.GGIW = GGIW_hat;
        
        
        tracks{i}(2:end) = [];
        
        
        
        
    end
end

table = ones(1,n_tt);
wAssoc = 0;

end