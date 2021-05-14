% Code by Yuxuan Xia

function [PPP,MBM] = predictPMBM(PPP,MBM,model)

% Predict existing PPP (for each GGIW component)
PPP.w = PPP.w + log(model.Ps);
% No trajectory information is kept for PPP (only for this implementation)
PPP.GGIW = arrayfun(@(x) predictGGIWPPP(x,model), PPP.GGIW);

% Incorporate PPP birth
PPP.w = [PPP.w;log(model.birth.w)];
PPP.GGIW = [PPP.GGIW;model.birth.GGIW];

% Predict MBM
n_track = length(MBM.track);
for i = 1:n_track % for each hypothesis tree (track)
    nh = length(MBM.track{i});
    for h = 1:nh % for each local hypothesis (Bernoulli)
        MBM.track{i}(h).Bern.GGIW = predictGGIW(MBM.track{i}(h).Bern.GGIW,model);
        MBM.track{i}(h).Bern.r = MBM.track{i}(h).Bern.r*model.Ps;
    end
end

end

