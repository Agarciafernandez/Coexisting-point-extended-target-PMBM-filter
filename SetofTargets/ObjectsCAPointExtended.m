% Code by Yuxuan Xia

function [partitions,lik,Nj,Berns,Liks] = ObjectsCAPointExtended(MBM,PPP,W,gating_matrix_d,gating_matrix_u,gating_matrix_up,model)

%number of measurements
m = size(W,2);

%generate multiple measurement partitions using DBSCAN
P = genPartitions_dbscan(W,model);

%number of partitions
np = length(P);

%convert measurement indices to boolean representation
C = cell(0,1);
nc = zeros(np,1);
for i = 1:np
    %number of measurement cells in the i-th measurement partition
    nc(i) = length(P{i});
    for j = 1:nc(i)
        temp = P{i}{j};
        P{i}{j} = false(m,1);
        P{i}{j}(temp) = true;
    end
    C = [C;P{i}'];
end

%find all the unique measurement cells
[unique_C,~,IC] = unique(cell2mat(C')','rows','stable');

%number of unique measurement cells
nuc = size(unique_C,1);

%for each measurement partition, find indices of measurement cell it contains
iP = cell(np,1);
idx = 0;
for i = 1:np
    iP{i} = sort(IC(idx+1:idx+nc(i)));
    idx = idx + nc(i);
end

%loop through each track (local hypothesis tree)

%number of pre-existing tracks
n_tt = length(MBM.track);
Cij = cell(n_tt,1);
Lij_meas = cell(n_tt,1);    %measurement update likelihood
Lij_miss = cell(n_tt,1);    %misdetection likelihood
for i = 1:n_tt
    %number of hypotheses in track i
    num_hypo = length(MBM.track{i});
    Cij{i} = cell(num_hypo,1);
    Lij_meas{i} = cell(num_hypo,1);
    Lij_miss{i} = zeros(num_hypo,1);
    
    Bern_meas{i} = cell(num_hypo,1);
    Bern_miss{i} = cell(num_hypo,1);
    
    for j = 1:num_hypo
        
        ii = 0;
        
        %for each local hypothesis, find indices of measurements inside its gate
        m_idx_ingate = gating_matrix_d{i}(:,j);
        %find measurement cells can be used to update this local hypothesis
        for k = 1:nuc
            nuck = sum(unique_C(k,:));
            if nuck == sum(unique_C(k,:).*m_idx_ingate')
                
                ii = ii + 1;
                
                Cij{i}{j} = [Cij{i}{j};k];
                %compute the association likelihood of assigning the kth
                %unique measurement cell to the jth local hypothesis in the
                %ith track
                [Bern_meas{i}{j}{ii},meas_lik] = detectionBernPointExtended(MBM.track{i}(j).Bern,W(:,unique_C(k,:)),model);
                Lij_meas{i}{j} = [Lij_meas{i}{j};meas_lik];
            end
        end
        [Bern_miss{i}{j},Lij_miss{i}(j)] = misdetectionBernPointExtended(MBM.track{i}(j).Bern,model);
    end
end

%compute the likelihood of being new tracks
n_ppp = length(PPP.w);
n_ppp_p= length(PPP.wp);
Lij_new = zeros(nuc,1);

Bern_new = cell(nuc,1);

for i = 1:nuc
    %for each measurement cell, determine whether it is inside the gate of
    %each PPP GGIW component
    in_gate = false(n_ppp,1);
    nuci = sum(unique_C(i,:));
    for j = 1:n_ppp
        if nuci == sum(unique_C(i,:).*gating_matrix_u(:,j)')
            in_gate(j) = true;
        end
    end
    if(size(W(:,unique_C(i,:)),2)>1)
        %It is an extended target
        if sum(in_gate) > 0
            [Bern_new{i},Lij_new(i)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,unique_C(i,:)),model);
            
        else
            %if not, treat all as clutter
            Lij_new(i) = sum(unique_C(i,:))*log(model.lambda_fa);
            Bern_new{i}.r = 0;
        end
    else
        %It can be a point target as well
        in_gate_point=false(n_ppp_p,1);
        for j = 1:n_ppp_p
            if nuci == sum(unique_C(i,:).*gating_matrix_up(:,j)')
                in_gate_point(j) = true;
            end
        end
        
        if sum(in_gate) > 0 || sum(in_gate_point) > 0
            [Bern_new{i},Lij_new(i)] = detectionPPP_joint(PPP.w(in_gate),PPP.GGIW(in_gate),PPP.wp(in_gate_point),PPP.Point(in_gate_point),...
                W(:,unique_C(i,:)),model);
        else
            %if not, treat all as clutter
            Lij_new(i) = sum(unique_C(i,:))*log(model.lambda_fa);
            Bern_new{i}.r = 0;
        end
    end
end

%loop through each global hypothesis
J = length(MBM.w);
Nj = zeros(J,1);
lik = [];
partitions = cell(J,1);
Berns = cell(J,1);
Liks = cell(J,1);
for j = 1:J
    %indices of tracks included in this global hypothesis
    track_indices = find(MBM.table(j,:)>0);
    %number of local hypotheses in this global hypothesis
    nj = length(track_indices);
    %for each local hypothesis, find indices of measurement cells that can
    %be associated to it and their corresponding association likelihoods
    cij = cell(nj,1);
    lij = cell(nj,1);
    lij_miss = zeros(nj,1);
    
    bernij = cell(nj,1);
    bernij_miss = cell(nj,1);
    
    for i = 1:nj
        cij{i} = Cij{track_indices(i)}{MBM.table(j,track_indices(i))};
        lij{i} = Lij_meas{track_indices(i)}{MBM.table(j,track_indices(i))};
        lij_miss(i) = Lij_miss{track_indices(i)}(MBM.table(j,track_indices(i)));
        
        bernij{i} = Bern_meas{track_indices(i)}{MBM.table(j,track_indices(i))};
        bernij_miss{i} = Bern_miss{track_indices(i)}{MBM.table(j,track_indices(i))};
    end
    %loop through each measurement partition
    for i = 1:np
        %construct cost matrix called by Murty
        C = Inf(nc(i),nj+nc(i));
        for i1 = 1:nc(i)
            for i2 = 1:nj
                %check if the i1-th measurement cell in the p-th partition
                %can be associated to the selected local hypothesis
                
                idx = ~logical(cij{i2} - iP{i}(i1));
                if any(idx)
                    C(i1,i2) = -(lij{i2}(idx) - lij_miss(i2));
                end
                
%                 idx = find(cij{i2} == iP{i}(i1),1);
%                 %if not empty, find its corresponding association weight
%                 if ~isempty(idx)
%                     %cost matrix for existing objects
%                     C(i1,i2) = -(lij{i2}(idx) - lij_miss(i2));
%                 end
            end
            
            %cost matrix for new objects
            C(i1,nj+i1) = -Lij_new(iP{i}(i1));
        end
        
        %obtain M best assignments using Murty's algorithm
        
        %Crouse's Murty
        %[col4rowBest,~,gainBest] = kBest2DAssign(C,ceil(exp(MBM.w(j)+log(model.M))));
        
        %Angel's Murty
        [col4rowBest,gainBest]= murty(C,ceil(exp(MBM.w(j)+log(model.M))));
        col4rowBest = col4rowBest';
        gainBest = gainBest';

        
        
        
        
        %restore weights
        lik = [lik;-gainBest + sum(lij_miss) + MBM.w(j)];
        %number of assignments
        Mh = length(gainBest);
        partition = cell(Mh,1);
        berns = cell(Mh,1);
        liks = cell(Mh,1);
        for h = 1:Mh
            partition{h} = cell(nj+nc(i),1);
            berns{h} = cell(nj+nc(i),1);
            liks{h} = cell(nj+nc(i),1);
            
            for ii = 1:nj
                berns{h}{ii} = bernij_miss{ii};
                liks{h}{ii} = lij_miss(ii);
            end
            
            for p = 1:nc(i)
                partition{h}{col4rowBest(p,h)} = find(unique_C(iP{i}(p),:)==true);
                if col4rowBest(p,h) <= nj
                    idx = find(cij{col4rowBest(p,h)} == iP{i}(p),1);
                    berns{h}{col4rowBest(p,h)} = bernij{col4rowBest(p,h)}{idx};
                    liks{h}{col4rowBest(p,h)} = lij{col4rowBest(p,h)}(idx);
                else
                    berns{h}{col4rowBest(p,h)} = Bern_new{iP{i}(p)};
                    liks{h}{col4rowBest(p,h)} = Lij_new(iP{i}(p));
                end
            end
            idx_empty = cellfun('isempty',partition{h});
            idx_empty(1:nj) = false;
            partition{h} = partition{h}(~idx_empty);
            berns{h} = berns{h}(~idx_empty);
            liks{h} = liks{h}(~idx_empty);
        end
        partitions{j} = [partitions{j};partition];
        Berns{j} = [Berns{j};berns];
        Liks{j} = [Liks{j};liks];
    end
    %number of new global hypotheses created by the j-th predicted global hypothesis
    Nj(j) = length(partitions{j});
end


