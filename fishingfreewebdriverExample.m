%% fishingfreewebdriverLocal.m

%% Author --------------------------------------------------------------
% name: Valentin Cocco
% mail: valentin.cocco@ens.fr
% creation: 4-2-2018

%% Description ---------------------------------------------------------
% Using the 1100-foodweb-structure stock, run simulations for 4000 timesteps without fishing effort locally.
% This serves as an example of the ecological dynamics of the aquatic
% networks without fishing. 

% Calls:
%   - webproperties.m
%   - setup_default.m
%   - differential.m

%% Updates -----------------------------------------------------------------
% when: 8-14-2019
% who: Paul Glaum (prglaum@umich.edu)
% what: used to run simulations without fishing effort for 500 different foodweb structures.

%% 1. UPLOADING OF THE NETWORK STRUCTURES

%Change to Data folder to load the Food-Web networks 'Webs1100.mat' file
cd('Data')
load('Webs1100.mat')
cd('..')

nSim=1; %number of simulations
spe=30; %number of species
t1=4000; %last timestep
tspan=0:t1;


    k=randi(500);    
    sprintf('Simulation %d/%d',k,nSim)
    
    web=webs{k,1};
    fish=webs{k,2};
    B0=webs{k,3};
    
    %% 2. CALCULATION OF THE INITIAL STRUCTURAL PROPERTIES
    [tmp, T]=webproperties(web);
    Fish=nnz(fish)/spe;
    properties=[tmp,Fish];
    
    %% 3. SET THE BIOLOGICAL PARAMETERS
    [r,K,y,e,Bs,c,h,ax_ar,Z,po,Bext]=setup_default(web,fish);
    x=ax_ar.*(Z.^(-po.*(T-1)));
    
    %% 4. ECONOMIC PARAMETERS
    mu=0;
    co=1;
    ca=0.01;
    a=1;
    b=0;
    price='linear'; 

    %% 5. ATN MODEL
    E0=[];
    harv=false(spe,1);
    X0=[B0,E0];
    options=odeset('RelTol',10^-8,'AbsTol',10^-8);
    [t,X] = ode45(@(t,X) differential(t,X,x,y,r,K,e,h,c,Bs,web,harv,mu,co,ca,a,b,price,Bext),tspan,X0,options);
    B=X(:,1:spe);
    E=X(:,spe+1:end);
    B(B<Bext)=0;
    E(E<0)=0;
    X=[B,E];
    

    %% 6. Plot Results
figure
set(gcf,'color','w');
%Time series of all trophic species
subplot(2,2,1)
plot(X(:,:));
%Time series of all trophic species that survive to end of simulation
subplot(2,2,2)
plot(X(:,~ext));
%Time series of harvested species, empty in this fishing free initial
%simulation
subplot(2,2,3)
plot(X(:,harv))
%Time series of all labeled fish species
subplot(2,2,4)
plot(X(:,fish))
