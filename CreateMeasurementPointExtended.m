function Z=CreateMeasurementPointExtended(model,object_tracks,K)


range_c =model.range_c;
Z = cell(K,1);
R=model.measmodel.R;  %R for point targets
for k = 1:K
    
    
    N_k=length(object_tracks);
    
    for i = 1:N_k
        
        t_birth=object_tracks(i).birthTime;
        t_death=object_tracks(i).deathTime;
        
        if(and(k>=t_birth,k<=t_death))
            %Target is alive
            if(isempty(object_tracks(i).g))
                %We have a point target
                x=object_tracks(i).x(1:2,k-t_birth+1);
                if rand < model.Pd              % simulate missed detection
                    p = mvnrnd(x,R,1);
                    Z{k} = [Z{k} p'];
                end
                
            else
                %We have an extended target
                x=object_tracks(i).x(1:2,k-t_birth+1);
                X=object_tracks(i).X(:,:,k-t_birth+1);
                
                if rand < model.Pd              % simulate missed detection
                    Nd=poissrnd(object_tracks(i).g);
                    p = mvnrnd(x,X,Nd); %No R in the measurement model
                    Z{k} = [Z{k} p'];
                end
            end
        end
    end
    % add false alarm
    N_c = poissrnd(model.lambda_c);
    C = repmat(range_c(:,1),[1 N_c])+ diag(range_c*[-1;1])*rand(2,N_c);
    Z{k} = [Z{k} C];
end