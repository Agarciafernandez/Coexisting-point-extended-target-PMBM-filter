% This code implements the multi-Bernoulli mixture (MBM) filter
% for coexisting point and extended targets (PMBM filter with Poisson
% intensity equal to zero, and multi-Bernoulli birth)

% A. F. Garcia-Fernandez, J. L. Williams, L. Svensson and Y. Xia, 
% "A Poisson multi-Bernoulli mixture filter for coexisting point and extended targets," 
% in IEEE Transactions on Signal Processing, 2021 doi: 10.1109/TSP.2021.3072006.

% The code is based on the code for the extended target
% PMBM filter developed by Yuxuan Xia
% https://github.com/yuhsuansia/Extended-target-PMBM-and-PMB-filters

% Authors: Angel F. Garcia-Fernandez and Yuxuan Xia


clear;
close all;
dbstop if error

addpath('Third-party code','Common','GOSPA extended target','SetofTargets')


%Choose data assocation method (only one available in this implementation): dataAssocMethod 1: conventional two-step
%approach clustering (DBSCAN) + assignment (MURTY)
dataAssocMethod = 1;
if dataAssocMethod == 1
    model.dataAssocMethod = 1;
else
    error('Please choose a valid data association method!')
end


%Choose a scenario: Scenario 1: 27 targets born at four different
%locations; Scenario 2: targets move in proximity (a broad birth prior).

%Scenario of the simulation
Scenario_point_extended_close;




%We consider MB birth
%If we only have one component for point and one component for extended
%we can do
model.birth.Bern.r=model.birth.w+model.birth.wp;
model.birth.Bern.GGIW=model.birth.GGIW(1);
model.birth.Bern.Point=model.birth.Point(1);
model.birth.Bern.Point.c=model.birth.wp/model.birth.Bern.r;

%We also remove the PPP birth component
model.birth.wp=[];
model.birth.w=[];
model.birth.GGIW=[];
model.birth.Point=[];





%True: (Track-oriented) MB implementation
%False: (Track-oriented) MBM implementation
model.ifTOPMB=false;



%If plot
ifplot = false;




%Number of Monte Carlo runs
Nmc=100;


%Number of time steps
K = model.K;
Nsteps=K;
%GOSPA errors across time
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error



%Initialise memory
estimates = cell(K,1);




%We create the groundTruth (It is just a reshaping of the information)
groundTruth=cell(K,1);
N_tracks=length(object_tracks);

for k=1:K
    N_k=0;
    for i = 1:N_tracks
        t_birth=object_tracks(i).birthTime;
        t_death=object_tracks(i).deathTime;
        if(and(k>=t_birth,k<=t_death))
            %Target is alive (objects are still alive at t_death)
            N_k=N_k+1;
            groundTruth{k}(N_k,1).x=object_tracks(i,1).x(:,k-t_birth+1);
            if(~isempty(object_tracks(i,1).X))
                %It is an extended target
                groundTruth{k}(N_k,1).X=object_tracks(i,1).X(:,:,k-t_birth+1);
            else
                %It is a point target
                groundTruth{k}(N_k,1).X=[];
            end
        end
    end
end



rand('seed',9)
randn('seed',9)


for i=1:Nmc
    
    %PPP initialisation
    PPP.w = log(model.birth.w); %weights in logarithm
    PPP.GGIW = [];  %We do not add these structures for MB birth
    PPP.wp = log(model.birth.wp);
    PPP.Point=[];
    
    %MBM initialisation
    MBM.w = [];     % Global hypotheses weights
    MBM.track = {}; % Local hypotheses trees
    MBM.table = []; % Global hypotheses look-up table
    
    
    MBM.w=0;
    MBM.table=1;
    MBM.track{1}(1,1).Bern=model.birth.Bern;
    MBM.track{1}(1,1).lik=0;
    MBM.track{1}(1,1).assocHistory(1,1).t=1;
    MBM.track{1}(1,1).assocHistory(1,1).meas=0;
    
    
    %Simulate measurements
    Z=CreateMeasurementPointExtended(model,object_tracks,K);
    
    
    %fprintf('Time step: ');
    
    tic
    for k = 1:K
        
        %Print info
        %fprintf('%d ',k);
        
        %Update step
        [PPPu,MBMu] = updatePMBMPointExtended(PPP,MBM,Z{k},k,model);
        
        
        
        
        %estimate of the multi-target states
        estimates{k} = estimatorPointExtended(MBMu,model);
        
        
        %We calculate the GOSPA error
        [gospa,~, decomp_cost] = GOSPA_point_extended(groundTruth{k},estimates{k},p,c,2);
        gospa_loc=decomp_cost.localisation;
        gospa_mis=decomp_cost.missed;
        gospa_fal=decomp_cost.false;
        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+gospa.^2;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_false_t_tot(k)=squared_gospa_false_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        
        
        %%%%%%
        if ifplot
            % For illustration purposes
            figure(2)
            clf
            plot(Z{k}(1,:),Z{k}(2,:),'kx','linewidth',1)
            hold on
            N_k=length(object_tracks);
            
            for j = 1:N_k
                
                t_birth=object_tracks(j).birthTime;
                t_death=object_tracks(j).deathTime;
                
                if(and(k>=t_birth,k<t_death))
                    %Target is alive
                    if(isempty(object_tracks(j).g))
                        %We have a point target
                        x=object_tracks(j).x(1:2,k-t_birth+1);
                        plot(x(1),x(2),'xg','linewidth',2)
                        
                    else
                        %We have an extended target
                        x=object_tracks(j).x(1:2,k-t_birth+1);
                        X=object_tracks(j).X(:,:,k-t_birth+1);
                        %plot(x(1),x(2),'b')
                        [cx,cy]=Sigmacircle(x(1),x(2),X,3);
                        h2 = plot(cx,cy,'b-','linewidth',2);
                        
                    end
                end
            end
            
            for j = 1:length(estimates{k})
                x=estimates{k}(j).x;
                if(isempty(estimates{k}(j).X))
                    %This is a point target
                    plot(x(1),x(2),'-xm','linewidth',2)
                    
                else
                    %this is an extended target
                    X=estimates{k}(j).X;
                    [cx,cy]=Sigmacircle(x(1),x(2),X,3);
                    h2 = plot(cx,cy,'r--','linewidth',2);
                end
                
            end
            xlim([-500,500])
            ylim([-500,500])
            grid on
            xlabel('x axis (m)')
            ylabel('y axis (m)')
            
        end
        
        PPP=PPPu;
        MBM=MBMu;
        
        %Prediction Step
        [PPP,MBM] = predictPMBMPointExtended_MB(PPP,MBM,model);
        
    end
    
    simulation_time = toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(simulation_time), ' sec'])
    
end


%Root mean square GOSPA errors at each time step
rms_gospa_t=sqrt(squared_gospa_t_tot/Nmc);
rms_gospa_loc_t=sqrt(squared_gospa_loc_t_tot/Nmc);
rms_gospa_false_t=sqrt(squared_gospa_false_t_tot/Nmc);
rms_gospa_mis_t=sqrt(squared_gospa_mis_t_tot/Nmc);


%Root mean square GOSPA errors across all time steps
rms_gospa_tot=sqrt(sum(squared_gospa_t_tot)/(Nmc*Nsteps))
rms_gospa_loc_tot=sqrt(sum(squared_gospa_loc_t_tot)/(Nmc*Nsteps))
rms_gospa_false_tot=sqrt(sum(squared_gospa_false_t_tot)/(Nmc*Nsteps))
rms_gospa_mis_tot=sqrt(sum(squared_gospa_mis_t_tot)/(Nmc*Nsteps))

%Save results
if(model.ifTOPMB)
    %This is PMB filter
    save(['MB_point_extended_pd_',int2str(100*model.Pd),'_clut_',int2str(lambda_c)])
else
    %This is PMBM filter
    save(['MBM_point_extended_pd_',int2str(100*model.Pd),'_clut_',int2str(lambda_c)])
    
    
end

figure(2)
plot(1:Nsteps,rms_gospa_t,'blue',...
    1:Nsteps,rms_gospa_loc_t,'r',...
    1:Nsteps,rms_gospa_false_t,'black',...
    1:Nsteps,rms_gospa_mis_t,'green','Linewidth',1.3)
grid on
xlabel('Time step')
ylabel('RMS GOSPA error')
legend('Total','Localisation','False','Missed')

