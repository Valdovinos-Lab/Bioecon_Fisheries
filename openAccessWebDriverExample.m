%LinearWebDriverExample.m

%% Author --------------------------------------------------------------
% name: Paul Glaum
% mail: prglaum@umich.edu
% creation: 7-10-2018

%% Description ---------------------------------------------------------
% 1. Run the second stage of simulation with Open Access fishing effort on the conserved foodwebs.

% Calls:
%   - webproperties.m
%   - setup_default.m
%   - differential.m
%   - save_data4.m 
% LinearWebDriverExample.m is set to use the linear pricing model described
% in the manuscript. For interrested researchers, there are different 
% options for this model. Options are: isoelastic (replace 'linear' with
% 'isoelastic') and non-linear, isoelastic (replace 'linear' with 'nl-ni').

%% LOAD DATA -------------------------------------------------------------------------------------
%This loads a file that contains the conserved webs from the fishing free films. 
cd('Data')
load('SimConsLin.mat')
cd('..')

%% SIMULATIONS -------------------------------------------------------------------------
tspan=0:4000;

    %% SETUP
    k=randi(length(SimConsLin));
    web=SimConsLin(k).topo.web; %randomly choose a conserved web from the list used in this study
    spe=length(web); % # of species? 
    fish=SimConsLin(k).topo.fish;  %delineate the fish
    ext=SimConsLin(k).free.ext; % delinate the fish that went extinct 
    B0=SimConsLin(k).free.B; %biomass after the free stage (first 4000 time steps) 
    
    T=SimConsLin(k).initial.Troph; %array of Trophic levels 
    [r,K,y,e,Bs,c,h,ax_ar,Z,po,Bext]=setup_default(web,fish);
    x=ax_ar.*(Z.^(-po.*(T-1))); %use the T to set up the metabolic rates
    mu=0.3; %mu: market reactivity
    co=.1; %co: cost of fishing per unit of effort
    ca=0.01; %ca: catchability per unit of effort per unit of biomass
    a=30;
    b=0.1;%0.01;
    price='linear';
   
    nharvest=1; % # of fish harvested
    harv=false(spe,1); %vector of 0s
    tmp=false(nnz(fish'==1 & ext==0),1); %nnz=# of non0 elements, why do we need the apostrophy? 
            %vector of 0s the length of the # of nonextinct fish
    
    %Choose either a random fish, or...
    %ind=randperm(nnz(fish'==1 & ext==0),nharvest); %randomly select one (non extinct) fish to harvest, index it
    %The fish with the highest population biomass
    [M,ind]=max(B0(fish'==1 & ext==0));
    
    tmp(ind)=true; %use that index to mark that spot in the 0s vector
    harv(fish'==1 & ext==0)=tmp; %assign that to the harv list
    
    %% RUN
    E0=1;%set E0 for the simulation starting point
    sprintf('Simulation %d/%d mu %d',k,length(SimConsLin),mu)
    X0=[B0,E0];
    options=odeset('RelTol',10^-8,'AbsTol',10^-8);
    [t,X] = ode45(@(t,X) differential(t,X,x,y,r,K,e,h,c,Bs,web,harv,mu,co,ca,a,b,price,Bext),tspan,X0,options);
    B=X(:,1:spe);
    E=X(:,spe+1:end);
    B(B<Bext)=0;
    E(E<0)=0;
    X=[B,E];
        
figure
set(gcf,'color','w');
%Time series of all trophic species
subplot(2,2,1)
plot(X(:,1:30));
%Time series of fishing Effort
subplot(2,2,2)
plot(X(:,31));
%Time series of harvested species
subplot(2,2,3)
plot(X(:,harv))
%Time series of all labeled fish species
subplot(2,2,4)
plot(X(:,fish))
