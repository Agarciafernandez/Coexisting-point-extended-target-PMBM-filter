% Code by Yuxuan Xia

function [Bern,lik] = misdetectionBernPointExtended(Bern,model)

% For extended target tracking, misdetection can either due to the fact
% that the target is not detected or the fact that the target is detected
% but generates 0 measurements. We should take both alternatives into
% account.

%Extended target update

%Probability that the target is alive and detected but generates 0 measurement
temp = model.Pd*(Bern.GGIW.b/(Bern.GGIW.b+1))^Bern.GGIW.a;
%Porbability that the target is alive but not detected
% model.Qd;

%Total misdetection probability
qD = model.Qd + temp;

%Normalised weights of different events that cause misdetection
w1 = model.Qd/qD;
w2 = temp/qD;

%No change for trajectory is not detected
GGIW1 = Bern.GGIW;

%Update Gamma parameters for target generates 0 measurement
GGIW2 = GGIW1;
GGIW2.b = GGIW2.b + 1;

%Misdetection likelihood
lik_extended = 1 - Bern.r + Bern.r*qD;

%Merge these two misdetection hypotheses
[~,Bern.GGIW.a,Bern.GGIW.b] = gammaMerge([w1;w2],[GGIW1.a;GGIW2.a],[GGIW1.b;GGIW2.b]);

%Point target part
if(isfield(Bern,'Point'))
    %There is a point target part
    c=Bern.Point.c;
    lik_total_prev=c*model.Qd+(1-c)*qD;
    lik= 1-Bern.r + Bern.r*lik_total_prev;
    
    Bern.r=Bern.r*lik_total_prev/lik;
    
    Bern.Point.c=c*model.Qd/lik_total_prev;
    
    lik=log(lik);
    
    
else
    %The target is extended. No point target component.
    %Updated Bernoulli existence probability
    Bern.r = Bern.r*qD/lik_extended;
    
    lik = log(lik_extended);
    
end
