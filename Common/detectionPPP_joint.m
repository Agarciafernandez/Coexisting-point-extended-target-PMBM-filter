% Code by Yuxuan Xia

function [Bern,lik] = detectionPPP_joint(we,GGIW,wp,Point,C,model)

%C only has one measurement so the new target may be point, extended or
%clutter

%Extended target update

%GGIW state update given measurement set
[GGIW_c,lik_c] = updateGGIWforPPP(GGIW,C,model);

%Weight of new local hypothesis
w_c = lik_c + we + log(model.Pd);

%Merge different components created by different PPP components
[w_hat,Bern.GGIW] = GGIW_merge_wrap(w_c,GGIW_c);

%Point target update

%GGIW state update given measurement set
[Point_c,lik_c] = updatePointforPPP(Point,C,model);

%Weight of new local hypothesis
w_c = lik_c + wp + log(model.Pd);

%Merging
w_c_exp=exp(w_c);
w_c_exp_norm=w_c_exp/sum(w_c_exp);
x_u=zeros(4,1);
P_u=zeros(4);

for i=1:length(Point_c)
    x_i=Point_c(1).xp;
    P_i=Point_c(1).Pp;
    
    x_u=x_u+w_c_exp_norm(i)*x_i;
    P_u=P_u+w_c_exp_norm(i)*(x_i*x_i'+P_i);
end

P_u=P_u-x_u*x_u';

lik_no_log=w_hat+sum(w_c_exp)+model.lambda_fa;
lik=log(lik_no_log);

Bern.r=(w_hat+sum(w_c_exp))/lik_no_log;
Bern.Point.xp=x_u;
Bern.Point.Pp=P_u;
Bern.Point.c=sum(w_c_exp)/(w_hat+sum(w_c_exp)); %Probability that the target is a point target
 


end

