% Code by Yuxuan Xia

function estimates = estimatorPointExtended(MBM,model)

%Number of estimated targets
N_k=0;

%Target extent dimension
d = 2;
estimates=cell(0);
%Extract estimates from the highest weight global hypothesis (multi-Bernoulli)
[~,idx] = max(MBM.w);
table_entry = MBM.table(idx,:);
for i = 1:length(table_entry)
    %Only extract estimates from Bernoullis with existence probability
    %larger than a pre-defined threshold
    if table_entry(i) > 0 && MBM.track{i}(table_entry(i)).Bern.r >= model.exist_r
        N_k=N_k+1;
        %Class Probability
        if(~isfield(MBM.track{i}(table_entry(i)).Bern,'Point'))
            %It is an extended target
            c=0;
        else
            
            c=MBM.track{i}(table_entry(i)).Bern.Point.c;
        end
        
        if(c>0.5)
            %Target estimate is a point target
            estimates(N_k,1).x=MBM.track{i}(table_entry(i)).Bern.Point.xp;
            estimates(N_k,1).X=[];
        else
            %Target estimate is an extended target
            
            estimates(N_k,1).x=MBM.track{i}(table_entry(i)).Bern.GGIW.m;
            estimates(N_k,1).X=MBM.track{i}(table_entry(i)).Bern.GGIW.V/(MBM.track{i}(table_entry(i)).Bern.GGIW.v-2*d-2);
            
        end
        %         estimates.g = [estimates.g MBM.track{i}(table_entry(i)).Bern.GGIW.a/MBM.track{i}(table_entry(i)).Bern.GGIW.b];
        %         estimates.x = [estimates.x MBM.track{i}(table_entry(i)).Bern.GGIW.m];
        %         X = MBM.track{i}(table_entry(i)).Bern.GGIW.V/(MBM.track{i}(table_entry(i)).Bern.GGIW.v-2*d-2);
        %         estimates.X = cat(3,estimates.X,X);
    end
end




