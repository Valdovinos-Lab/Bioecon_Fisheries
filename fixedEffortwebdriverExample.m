%% constantwebdriverExample.m

%% Author --------------------------------------------------------------
% name: Valentin Cocco
% mail: valentin.cocco@ens.fr
% creation: 5-16-2018
% editor: Paul Glaum

%% Description ---------------------------------------------------------
% 1. Run the second stage of simulation with Fixed fishing effort on conserved foodwebs.

% Calls:
%   - webproperties.m
%   - setup_default.m
%   - differential.m

%% Updates -----------------------------------------------------------------
% when: 8-14-2019
% what: adaptation to updated analysis

%% LOAD DATA -------------------------------------------------------------------------------------
%This loads a file that contains the conserved webs from the fishing free films. 
cd('Data')
load('SimCons.mat')
cd('..')


%% SIMULATIONS -------------------------------------------------------------------------
tspan=0:4000;

    %% SETUP
    k=randi(length(SimCons)); 
    web=SimCons(k).topo.web; %randomly choose a conserved web from the list used in this study
    spe=length(web);
    fish=SimCons(k).topo.fish;
    ext=SimCons(k).free.ext;
    B0=SimCons(k).free.B;
    
    T=SimCons(k).initial.Troph;
    [r,K,y,e,Bs,c,h,ax_ar,Z,po,Bext]=setup_default(web,fish);
    x=ax_ar.*(Z.^(-po.*(T-1)));
    
    %mu=0, this keeps the Effort static at E0, the initial Effort value
    mu=0;
    
    co=1;
    ca=0.01;
    a=1;
    %b=0, this keeps the price static despite changing yield
    b=0;
    price='linear';
    nharvest=1;
    %randomly choose fish to harvest
    harv=false(spe,1);
    tmp=false(nnz(fish'==1 & ext==0),1); 
    ind=randperm(nnz(fish'==1 & ext==0),nharvest);
    tmp(ind)=true;%use that index to mark that spot in the 0s vector
    harv(fish'==1 & ext==0)=tmp;%assign that to the harv list
    
    
    %% EFFORT LEVELS
    Effort=[0,1,2,3,4,5,6,7,8,9,10,12,15,20]; % list of fixed effort levels
    
    %% RUN
    E0=datasample(Effort,1); %randomly choose a Fixed Effort level, change this to look at specific levels
    sprintf('Simulation %d/%d Effort %d',k,length(SimCons),E0)
    X0=[B0,E0];
    options=odeset('RelTol',10^-8,'AbsTol',10^-8);
    [t,X] = ode45(@(t,X) differential(t,X,x,y,r,K,e,h,c,Bs,web,harv,mu,co,ca,a,b,price,Bext),tspan,X0,options);
    B=X(:,1:spe);
    E=X(:,spe+1:end);
    B(B<Bext)=0;
    E(E<0)=0;
    X=[B,E];
        

%% PLOT RESULTS
figure
set(gcf,'color','w');
%Time series of all trophic species
subplot(2,2,1)
plot(X(:,1:30));
%Time series of Effort (unchanging in the Fixed Effort simulation)
subplot(2,2,2)
plot(X(:,31));
%Time series of harvested species
subplot(2,2,3)
plot(X(:,harv))
%Time series of all labeled fish species
subplot(2,2,4)
plot(X(:,fish))
