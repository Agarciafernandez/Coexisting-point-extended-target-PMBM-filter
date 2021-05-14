%Scenario for coexisting point and extended targets
rand('seed',10)
randn('seed',10)

K=100;
model.K = K;

% Effective window length for the gamma prediction
w_e_gamma = 20;
% Effective window length for the extent prediction
w_e_extent = 10;

% Parameters to control the agility of the prediction
model.tao = 1/(log(w_e_extent)-log(w_e_extent-1));
model.eta = 1/(1-1/w_e_gamma);

model.Ts = 1;   %sampling interval
sigma_v = 0.5;  %standard deviation of motion noise
sigma_r = 1;  %standard deviation of measurement noise


% Linear motion and measurement models
model.motionmodel = motionmodel2.cvmodel(model.Ts,sigma_v);
model.measmodel = measmodel.cvmeasmodel(sigma_r);


%Choose object survival probability
model.Ps = 0.99;


% Target detection probability
model.Pd = 0.95;




%Choose range of surveillance area
sensor_model.range_c = [-500 500;-500 500];
model.range_c = [-500 500;-500 500];


lambda_c=8; %False alarm rate (notice difference with lambda_fa)
model.lambda_c=lambda_c;
model.lambda_fa= lambda_c/prod(model.range_c(:,2)-model.range_c(:,1));

%Choose gamma distribution parameter
alpha = 40;
beta = 4;
sensor_model.gamma = makedist('Gamma','a',alpha,'b',1/beta);



%Specify inverse-Wishart parameters
model.motionmodel.iwish.V = 200*eye(2);
model.motionmodel.iwish.v = 20;

%Create ground truth model
nbirths = 1;
Pb=diag([200^2,200^2,4^2,4^2]);
birth_model = repmat(struct('w',log(0.06),'x',[],'P',Pb),[1,nbirths]);
birth_model(1).x  = [ 0; 0; 0; 0 ];

%We add the intensity for point targets

birth_model.wp=log(0.03);
birth_model(1).xp=[0;0;0;0];
birth_model(1).Pp=Pb;


model.birth.w = exp(birth_model.w);%weights of birth components

% Each GGIW component is described by the sufficient statistics of Gamma
% distribution (a,b), Gaussian distirbution (m,P) and inverse-Wishart
% distributino (v,V). For inverse-Wishart distribution, v is a scalar
% controls the distribution uncertainty. The larger v, the smaller
% uncertainty.
model.birth.GGIW = repmat(struct('a',5e3,'b',1e3,'m',[],'P',birth_model(1).Pp,'v',14,'V',10*eye(2)),[nbirths,1]);
for i = 1:nbirths
    model.birth.GGIW(i).m =birth_model(i).x;
end
%For Point targets
model.birth.wp=exp(birth_model.wp);
model.birth.Point(1).xp=birth_model(1).xp;
model.birth.Point(1).Pp=birth_model(1).Pp;

% Gating parameters
Pg = 0.999; %gating size in probability
Pg = 0.9999; %Change made by Angel

model.gamma= chi2inv(Pg,model.measmodel.d);
% Effective missed detection probability after applying gating
model.Qd = 1 - model.Pd*Pg;

% Pruning thresholds
model.threshold_r = 1e-3;   %existence probability of Bernoulli component
model.threshold_u = 1e-3;   %weight of mixture component in PPP
model.threshold_w = 1e-3;   %weight of global hypothesis (multi-Bernoulli)
%model.threshold_s = 1e-4;   %weight of trajectory is alive if exists

%model.recycle = 1e-1;       %recycling threshold
%model.merge = 4;            %merge threshold used to merge similar GGIWs
model.M = 20;              %cap of number of MBM components in PMBM
%model.num_iterations = 10;  %controls the number of iterations used in SO
%model.max_repetition = 2;   %controls the number of iterations used in SO

%extract target state from Bernoulli components with existence probability
%larger than this threshold
model.exist_r = 0.5;

% DBSCAN parameters, a grid search for hyperparameters
model.max_dist = 12;
model.min_dist = 0.1;
model.grid_dist = 0.1;

%Generate ground truth set of trajectories
[object_tracks] = trackgen(K,model.motionmodel,sensor_model,birth_model,model.Ps);

%Parameters used in GOSPA metric
c = 10;
p = 2;


%Plot the scenario
figure(1)
box on
grid on
hold on
cols = parula(length(object_tracks));
for it = 1:length(object_tracks)
    xx = object_tracks(it).x(1,:);
    yy = object_tracks(it).x(2,:);
    plot(xx,yy,'b','linewidth',1.3)
    plot(xx(1),yy(1),'xr','linewidth',2)
    text(xx(1)+15,yy(1),num2str(object_tracks(it).birthTime))
    plot(xx(11:10:end),yy(11:10:end),'xb','linewidth',2)

    if(~isempty(object_tracks(it).X))
        for ii = 1:10:size(xx,2)
            %illustrate the 3-sigma level of ellipse
            [cx,cy]=Sigmacircle(xx(ii),yy(ii),object_tracks(it).X(:,:,ii),3);
            %plot(cx,cy,'-','color',cols(it,:),'linewidth',2);
            plot(cx,cy,'-','color','b','linewidth',1.3);
        end
    end
end
xlabel('x axis (m)')
ylabel('y axis (m)')
grid on
axis([-400 400 -400 400])






function [object_tracks] = trackgen(K,motionmodel,sensormodel,birthmodel,P_S)

%surveillance area
range_c = sensormodel.range_c;


n = 0;

xmid=[range_c(1,2)+range_c(1,1);range_c(2,2)+range_c(2,1);0;0]/2;
Pmid = diag([36 36 16 16]);
kmid=fix(K/2);

F=motionmodel.F(xmid);
chol_Q=chol(motionmodel.Q)';

%Convert birthmodel.w
[log_w,log_sum_w] = normalizeLogWeights([birthmodel.w]);

[log_wp,log_sum_wp] = normalizeLogWeights([birthmodel.wp]);



for k = 1:kmid-10 %These are possible birth times
    %randomly sample number of births (extended targets)
    nb = poissrnd(exp(log_sum_w));
    
    if(k==1)
        nb=2;
    else
        nb=0;
    end
    
    
    for i = 1:nb
        %randomly sample birth component
        idx = find(rand<cumsum(exp(log_w)),1,'first');
        n = n+1;
        %randomly sample kinematic state
        
        
        
        
        
        object_tracks(n,1).x = mvnrnd(xmid, Pmid)';
        %randomly sample Poisson rate
        object_tracks(n,1).g = random(sensormodel.gamma);
        %randomly sample extent matrix
        %Note that for coordinated turn motion model, it is common to
        %rotate the extent matrix using the rotation matrices R obtained
        %from the object heading, such that X = R'XR
        Matrix_extent=iwishrnd(motionmodel.iwish.V, motionmodel.iwish.v);
        object_tracks(n,1).X = Matrix_extent;
        %birth time
        object_tracks(n,1).birthTime = k;
        %death time
        object_tracks(n,1).deathTime = kmid;
        
        time = kmid;
        xk = object_tracks(n,1).x;
        
        xpos = xmid(1);
        ypos = xmid(2);
        
        %We sample forward
        ifalive = 1;
        while ifalive && (xpos>=range_c(1,1))&&(xpos<=range_c(1,2))...
                &&(ypos>=range_c(2,1))&&(ypos<=range_c(2,2))&&(time+1<=K)
            
            %simulate the kinematic state
            xk = mvnrnd(motionmodel.f(xk), motionmodel.Q)';
            object_tracks(n,1).x = [object_tracks(n,1).x xk];
            object_tracks(n,1).deathTime = object_tracks(n,1).deathTime + 1;
            
            %We consider a constant extent matrix
            
            object_tracks(n,1).X = cat(3,object_tracks(n,1).X,...
                Matrix_extent);
            
            xpos = xk(1);
            ypos = xk(2);
            
            time=time+1;
            
            if rand > P_S
                ifalive = 0;
            end
        end
        
        %We now sample backwards
        xk=object_tracks(n,1).x(:,1);
        for kk=kmid-1:-1:k
            xk = F\(xk + chol_Q*randn(4,1));
            xpos = xk(1);
            ypos = xk(2);
            if (xpos<range_c(1,1))||(xpos>range_c(1,2))...
                    ||(ypos<range_c(2,1))||(ypos>range_c(2,2))
                object_tracks(n,1).birthTime=kk+1;          
                break;
            end
            
            
            object_tracks(n,1).x = [xk object_tracks(n,1).x];
            object_tracks(n,1).X = cat(3,Matrix_extent,object_tracks(n,1).X);
            
            
            
        end
        
    end
    
    %randomly sample number of births (point targets)
    
    nb = poissrnd(exp(log_sum_wp));
    
    if(k==5)
        nb=1;
    elseif(k==10)
        nb=1;
    else
        nb=0;
    end
    
    
    for i = 1:nb
        %randomly sample birth component
        idx = find(rand<cumsum(exp(log_wp)),1,'first');
        n = n+1;
        %randomly sample kinematic state
        object_tracks(n,1).x = mvnrnd(xmid, Pmid)';
        %No Poisson rate
        object_tracks(n,1).g = [];
        %No object extent
        object_tracks(n,1).X = [];
        %birth time
        object_tracks(n,1).birthTime = k;
        %death time
        object_tracks(n,1).deathTime = k;
        
        time = kmid;
        xk = object_tracks(n,1).x;
        
        xpos = xmid(1);
        ypos = xmid(2);
        
        ifalive = 1;
        while ifalive && (xpos>=range_c(1,1))&&(xpos<=range_c(1,2))...
                &&(ypos>=range_c(2,1))&&(ypos<=range_c(2,2))&&(time+1<=K)
            
            %simulate the kinematic state
            xk = mvnrnd(motionmodel.f(xk), motionmodel.Q)';
            object_tracks(n,1).x = [object_tracks(n,1).x xk];
            object_tracks(n,1).deathTime = object_tracks(n,1).deathTime + 1;
            
            xpos = xk(1);
            ypos = xk(2);
            
            time=time+1;
            
            if rand > P_S
                ifalive = 0;
            end
        end
        %We now sample backwards
        xk=object_tracks(n,1).x(:,1);
        for kk=kmid-1:-1:k
            xk = F\(xk + chol_Q*randn(4,1));
            xpos = xk(1);
            ypos = xk(2);
            if (xpos<range_c(1,1))||(xpos>range_c(1,2))...
                    ||(ypos<range_c(2,1))||(ypos>range_c(2,2))
                object_tracks(n,1).birthTime=kk+1;          
                break;
            end
            object_tracks(n,1).x = [xk object_tracks(n,1).x];
        end
        
        
    end
    
end

if n==0
    object_tracks = [];
end

end