% Code by Yuxuan Xia

function [Bern,lik] = detectionPPPPoint(wp,Point,C,model)

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

Bern.Point.xp=x_u;
Bern.Point.Pp=P_u;

Bern.r=sum(w_c_exp)/(sum(w_c_exp)+model.lambda_fa);
lik=log(sum(w_c_exp)+model.lambda_fa);

end

