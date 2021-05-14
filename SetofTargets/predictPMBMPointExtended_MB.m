% Code by Yuxuan Xia

function [PPP,MBM] = predictPMBMPointExtended_MB(PPP,MBM,model)

%This code considers MB birth

%PPP for extended targets
% Predict existing PPP (for each GGIW component)
PPP.w = PPP.w + log(model.Ps);
% No trajectory information is kept for PPP (only for this implementation)
PPP.GGIW = arrayfun(@(x) predictGGIWPPP(x,model), PPP.GGIW);

%PPP for point targets
PPP.wp=PPP.wp+log(model.Ps);


if(~isempty(PPP.Point))
F = model.motionmodel.F(PPP.Point(1).xp);
Q= model.motionmodel.Q;

for i=1:length(PPP.Point)
    PPP.Point(i).xp=F*PPP.Point(i).xp;
    PPP.Point(i).Pp=F*PPP.Point(i).Pp*F'+Q;
end
else
    %We need to obtain F and Q
    F = model.motionmodel.F( zeros(4,1)); %We just particularise for the 4 dimensional case.
Q= model.motionmodel.Q;
end
    


% Incorporate PPP birth
PPP.w = [PPP.w;log(model.birth.w)];
PPP.GGIW = [PPP.GGIW;model.birth.GGIW];
PPP.wp = [PPP.wp;log(model.birth.wp)];
PPP.Point=[PPP.Point;model.birth.Point];



% Predict MBM
n_track = length(MBM.track);
for i = 1:n_track % for each hypothesis tree (track)
    nh = length(MBM.track{i});
    for h = 1:nh % for each local hypothesis (Bernoulli)
        MBM.track{i}(h).Bern.GGIW = predictGGIW(MBM.track{i}(h).Bern.GGIW,model);
        MBM.track{i}(h).Bern.r = MBM.track{i}(h).Bern.r*model.Ps;
        
        if(isfield(MBM.track{i}(h).Bern,'Point'))
            x_ih= MBM.track{i}(h).Bern.Point.xp;
            P_ih= MBM.track{i}(h).Bern.Point.Pp;
            
            MBM.track{i}(h).Bern.Point.xp=F*x_ih;
            MBM.track{i}(h).Bern.Point.Pp=F*P_ih*F'+Q;
            %The class probability does not change, so we do not do anything     
        end
        
    end
end


%We add the new tracks
MBM.table=[MBM.table,ones(size(MBM.table,1),1)];

MBM.track{n_track+1}(1).Bern=model.birth.Bern;
MBM.track{n_track+1}(1).lik=0;
MBM.track{n_track+1}(1).assocHistory(1,1).t=1;
MBM.track{n_track+1}(1).assocHistory(1,1).meas=0;

end

