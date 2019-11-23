%% setup_default.m

%% Previous authors -------------------------------------------------------
% Program by: Rosalyn Rael
% Modified by Barbara Bauer (April 2009)

%% Last update ------------------------------------------------------------
% when: 4-2-2018
% who: Valentin Cocco
% mail: valentin.cocco@ens.fr
% what: personal writing for constant effort framework

%% Description ------------------------------------------------------------
% Set the biodynamical DEFAULT parameters of the foodweb.
% Used for  the constant effort framework
% Inputs:
%   - nicheweb: adjacency matrix (row i eats column j)
%   - fish: logical array; indicate which species are fishes
% Called by:
%   - sensitivity_check.m

%%
function [r,K,y,e,Bs,c,h,ax_ar,Z,po,Bext]=setup_default(nicheweb,fish)

%% ------------------------------------------------------------------------
% 1. DYNAMICAL PARAMETERS -------------------------------------------------
% -------------------------------------------------------------------------
% r: growth rate of i (non zero for basal species only)
% K: global carrying capacity of i (for basal species only)
% y: maximum consumption rate of i eating j
% e: assimilation efficiency of i eating j
% Bs: half saturation biomass for i eating j
% c: interference of predator i eating j with other predators of j
% h: Hill coefficient

%% r: growth rate ---------------------------------------------------------------------------
% constant r=1 for every basal species
spe=length(nicheweb);
r=zeros(spe,1);
basal=(sum(nicheweb,2)==0);
r(basal)=1;

%% K: global carrying capacity ----------------------------------------------------------------
K=ones(spe,1);

%% x: metabolic rate --------------------------------------------------------------------------
nonbasal=(sum(nicheweb,2)~=0);

ax_ar=zeros(spe,1);
ax_ar(nonbasal)=0.314;
ax_ar(fish)=0.88;

Z=100;

po=0.25;

%% y: maximal consumption rate per unit or metabolic rate -------------------------------------
y=zeros(spe,1);

y(nonbasal)=10;

%% e: assimilation efficiency -----------------------------------------------------------------
% carnivorous/herbivorous
%i eats non basal (carnivorous)
e=0.85*ones(spe);
%i eats basal (herbivorous)
e(:,basal)=0.45;

%% Bs: half saturation biomass ---------------------------------------------------------------
Bs=0.2*ones(spe,1);

%% c: feeding interference coefficient -------------------------------------------------------
% Holling-type II or III
c=zeros(spe,1);

%% h: Hill coefficient ------------------------------------------------------------------------
h=1.2;

%% Bext: extinction threshold -----------------------------------------------------------------
Bext=10^-6;
