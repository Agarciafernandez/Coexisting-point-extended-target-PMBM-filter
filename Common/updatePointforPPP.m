% Code by Yuxuan Xia

function [Point_c,lik_c] = updatePointforPPP(Pointp,z,model)




%Perform Kalman filter update for each mixture component

np = length(Pointp);
if(np==0)
    Point_c=[];
    lik_c=[];
    return;
    
else


H = model.measmodel.H(Pointp(1).xp);
R=model.measmodel.R;
lik_c = zeros(np,1);
Point_c = repmat(struct('xp',[],'Pp',[]),[np,1]);
Nz=2;

for i=1:np
    x_i= Pointp(1).xp;
    P_i=Pointp(1).Pp;
    
    %Kalman filter update
    S=H*P_i*H'+R;
    z_pred=H*x_i;
    K_gain=P_i*H'/S;
    x_u=x_i+K_gain*(z-z_pred);
    P_u=P_i-K_gain*H*P_i;

    %Likelihood goes in log
    Vs= chol(S);
    log_det_S= 2*log(prod(diag(Vs)));
    maha=(z-z_pred)'/S*(z-z_pred);
    
    
    lik_c(i)=-0.5*maha-1/2*log_det_S-Nz*log(2*pi)/2;
    %lik_c(i)=mvnpdf(z,z_pred,S);
    Point_c(i).xp=x_u;
    Point_c(i).Pp=P_u;
    
end
end

