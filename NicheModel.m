%% NicheModel.m

%% Previous authors -------------------------------------------------------
% Program by: Rosalyn Rosalyn Rael
% Barbara Bauer added a part removing species not connected to a basal 
% and changed the output to a 'row eats column' type matrix
% Coralie Picoche changed the 'removing part' in order not to change the connectance
% was adding the loop !
% Ref: Williams and Martinez, Nature, 2000.

%% Last update ------------------------------------------------------------
% when: 3-6-2018
% who: Valentin Cocco
% mail: valentin.cocco@ens.fr
% what: personal rewriting

%% Description ------------------------------------------------------------
% Produce a niche model foodweb
% Inputs:
%   - spe: number of species
%   - con: connectance (C=L/S^2 ; L: number of links)
%   - err: relative error on con
% Output:
%   - nicheweb: adjacency matrix (nicheweb(i,j)=1 <=> row i eats colomn j)
% Calls:
%   - isConnected.m
% Called by:
%   - webdriver.m

%% Function OR Script
function [nicheweb]= NicheModel(spe,con,err)
%spe=30; %(Berlow2009)
%con=0.15; %(Berlow2009)
%err=0.025; %(Williams2000: 0.03)

%% Safety conditions
tries=10000; %nb of tries to find a web which respect the fixed conditions
ok=0; %logical (0=false ; true=1) which indicates if the network is correct
% Check that the range of connectance include integer numbers of links
Lmin=ceil(con*(1-err)*spe^2);   %ceil : return the nearest superior integer
Lmax=floor(con*(1+err)*spe^2);  %floor : return the nearest inferior integer
if Lmin>Lmax
    error(['Impossible to create a foodweb with the given number of species, connectance +/- error. ', 'Please change the inputs.'])
end

%% Main loop
while (ok==0 && tries>0)
    ok=1;
    tries=tries-1;
    
    %% Assign a 1D-niche to each species
    n=rand(spe,1); %assign a niche value n to each species (uniform distribution)
    % Assign a range value r to each species (beta distribution)
    a=1;
    b=(1-2*con)/(2*con); %(E(r)=con)
    r=betarnd(a,b,spe,1);
    r=r.*n;
    % Assign the center of the niche c to each species (uniform distribution)
    c=min(rand(spe,1).*(n-r./2)+r./2 , 1-r./2);
    % Force the species with the smallest niche value to be basal (no prey)
    [~,idx] = sort(n); %ascending order
    r(idx(1))=0;
    
    %% Construct the adjacency matrix (if n(j) is included in the niche range of i, then i eats j)
    r_min=c-r./2;
    r_max=c+r./2;
    r_min=r_min*ones(1,spe);
    r_max=r_max*ones(1,spe);
    tmp=ones(spe,1)*n';
    nicheweb=(tmp>=r_min & tmp<=r_max); %/!\ a species i with r(i)=0 may eat j if c(i)=c(j)...
    
    %% Check the network is ok :
    %   - no isolated species
    %   - all species are (indirectly) connected to a basal
    %   - one unique network (no independant food web)
    %   - connectance = con +/- err
    %% 1. Is there any isolated species?
    predator=sum(nicheweb,1)'; %number of predator per species
    prey=sum(nicheweb,2); %number of prey per species
    link=predator+prey; %number of links per species
    %% 
    if all(link)==0
        ok=0;
    else
        %% 2. Is there any species not connected to a basal species ? 
        tot=zeros(spe,1); %tot(i)=1 : species i connected to a basal species
        basal=find(sum(nicheweb,2)==0);
        tmp=basal; %tmp: last species connected to basal species found
        tmpweb=nicheweb;
        tot(tmp)=1;
        while isempty(tmp)==0 %when tmp=[], there is no more species connected to basal ones to find
            [i,~]=find(tmpweb(:,tmp)~=0);
            tmpweb(:,tmp)=0;
            tmpweb(tmp,:)=0; %the last species "disappear" from the web to avoid infinite loops (i eats j, that eats k, that eats i...)
            tmp=unique(i); %i predators of tmp --> i are connected to basal species
            tot(tmp)=1;       
        end
        %% 
        if all(tot)==0 %at least one species is not connected to a basal
            ok=0;
        else
            %% 3. Is there more than one independant network ?
            z=isConnected(nicheweb);
            %% 
            if z==0
                ok=0;
            else
                %% 4. Is the actual connectance close enough from con ?
                con_web=nnz(nicheweb)/spe^2;
                con_min=(1-err)*con;
                con_max=(1+err)*con;
                %% 
                if con_web<con_min || con_web>con_max
                    ok=0;
                end
            end
        end
    end
end
%%
if tries==0
    error('Impossible to create a foodweb within 10,000 tries')
end