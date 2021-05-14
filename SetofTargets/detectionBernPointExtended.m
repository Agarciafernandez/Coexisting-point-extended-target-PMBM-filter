% Code by Yuxuan Xia

function [Bern,lik] = detectionBernPointExtended(Bern,C,model)



%Number of measurements
n_m=size(C,2);

%local hypothesis state update
[Bern.GGIW,lik_extended] = updateGGIW(Bern.GGIW,C,model);


if (~isfield(Bern,'Point')) 
    %local hypothesis weight
    lik = lik_extended + log(Bern.r) + log(model.Pd);
elseif (n_m>1 && isfield(Bern,'Point'))
     Bern = rmfield(Bern,'Point'); % We remove the field Point as it is an extended target
    %local hypothesis weight
    lik = lik_extended + log(Bern.r) + log(model.Pd);  
else
    z=C;
    %There is one measurement and we apply a KF update
    xp=Bern.Point.xp;
    Pp=Bern.Point.Pp;
    c=Bern.Point.c; %Point target probability
    H=[eye(2),zeros(2)];
    Nz=2;
    R=model.measmodel.R;
    
    S=H*Pp*H'+R;
    K=Pp*H'/S;
    z_pred=H*xp;
    
    xu=xp+K*(z-z_pred);
    Pu=Pp-K*H*Pp;
    
    Bern.Point.xp=xu;
    Bern.Point.Pp=Pu;
    
     %Likelihood goes in log
    Vs= chol(S);
    log_det_S= 2*log(prod(diag(Vs)));
    maha=(z-z_pred)'/S*(z-z_pred);
    
    
    lik_point=-0.5*maha-1/2*log_det_S-Nz*log(2*pi)/2;
    %lik_point=mvnpdf(z,z_pred,S);
    
    lik_exp=Bern.r*model.Pd*( exp(lik_point)*c+exp(lik_extended)*(1-c));
    
    %Updated class probability
    Bern.Point.c=Bern.r*model.Pd*exp(lik_point)*c/lik_exp;
    
    lik=log(lik_exp);
    

end


Bern.r = 1;


end
