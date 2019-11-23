%% webproperties.m

%% Previous authors
% Program by: Perrine Tonin and Coralie Pioche
% Ref:
%   - Dunne et al., 2006, PLOS Biolgy.
%   - Williams & Martinez, 2008, Journal of Animal Ecology.
%   - Williams & Martinez, 2000, Nature.

%% Last update
% when: 4-2-2018
% who: Valentin Cocco
% mail: valentin.cocco@ens.fr
% what: cannibalism included for CLust coefficient

%% Description
% Calculate the structural properties of a foodweb
% Inputs:
%   - nicheweb: matrix of adjacency. nicheweb(i,j)=1 <=> row i eats colomn j
% Outputs:
%   - properties: array containing 17 strucutural properties
%   - T: array containing trophic level of each species
% Called by:
%   - webdriver.m

%% 
function [properties, T]=webproperties(nicheweb)

%% 1. Bas: % of basal species (no prey)
spe=length(nicheweb);
basal=find(sum(nicheweb,2)==0);
Bas=length(basal)/spe;

%% 2. Top: % of top species (no predator)
top=find(sum(nicheweb,1)==0);
Top=length(top)/spe;

%% 3. Int: % of intermediate species
Int=1-Top-Bas;

%% 4. Can: % of cannibal species
%cannibal=diag(nicheweb)==1; %/!\ error with diag
cannibal=diag(nicheweb);
Can=nnz(cannibal)/spe;

%% 5. Herb: % of herbivore species (eat only basal species)
tmp1=sum(nicheweb(:,basal),2)~=0; %select species which eat basal species
nonbasal=find(sum(nicheweb,2)~=0);
tmp2=sum(nicheweb(:,nonbasal),2)==0; %select species which do not eat nonbasal species
herb= tmp1 & tmp2; %intersection
Herb=nnz(herb)/spe;

%% 6. MaxSim: mean of maximum similarity between species
% Remark: cannibalism included
S=zeros(spe); %matrix of similarity (symetric) : s_i_j: % of common links between i and j
for i=1:spe
    for j=1:spe
        if i~=j
            common_prey=nnz((nicheweb(i,:)==1 & nicheweb(j,:)==1));
            common_pred=nnz((nicheweb(:,i)==1 & nicheweb(:,j)==1));
            tot_prey=nnz((nicheweb(i,:)==1 | nicheweb(j,:)==1));
            tot_pred=nnz((nicheweb(:,i)==1 | nicheweb(:,j)==1));
            S(i,j)=(common_prey+common_pred)/(tot_prey+tot_pred);
        end
    end
end
Smax=max(S); %row vector containing the maximum element from each column
MaxSim=mean(Smax);

%% 7. VulSD: standard deviation of vulnerability (nb of predators)
vul=sum(nicheweb,1);
VulSD=std(vul);

%% 8. GenSD: standard deviation of generalism (nb of preys)
gen=sum(nicheweb,2);
GenSD=std(gen);

%% 9. LinkSD: standard deviation of degree (number of links)
LinkSD=std(vul'+gen);

%% 10. Clust: mean of % of realized links within the neighbours of k (k excluded)
% Watts & Strogatz 1998. /!\ here: directed networks.
clust=zeros(spe,1);

% a. without cannibalism
% We have to ignore neighbours whose k is their only neighbour (degree=1)
% nocanweb=nicheweb;
% nocanweb(logical(eye(spe)))=0;
% for k=1:spe
%     neighbours=find(nocanweb(k,:)'==1 | nocanweb(:,k)==1); %prey and predators of k (neighbours)
%     nk=length(neighbours);
%     tmpcluster=nocanweb(neighbours,neighbours); %adjacency matrix with k's neighbours only (k excluded)
%     clust(k)=nnz(tmpcluster)/(nk*(nk-1)); % #realized links/#potential links
% end

% b. with cannibalism
nocanweb=nicheweb;
nocanweb(logical(eye(spe)))=0;
for k=1:spe
    neighbours=find(nocanweb(k,:)'==1 | nocanweb(:,k)==1);
    nk=length(neighbours);
    tmpcluster=nicheweb(neighbours,neighbours);
    clust(k)=nnz(tmpcluster)/nk^2;
end

Clust=mean(clust, 'omitnan');

%% 11. Path: mean of the shortest path between each pair of species
% Watts & Strogatz 1998. /!\ here: undirected networks.
distweb=ones(spe); %(distance matrix betwenn nodes)
distweb(max(nicheweb,nicheweb')==0)=10^6; %(no link : infinite distance)
for k=1:spe
    for i=1:spe
        for j=1:spe
            distweb(i,j)=min(distweb(i,j),distweb(i,k)+distweb(k,j));
        end
    end
end
distweb=triu(distweb);
% distweb(logical(eye(spe)))=0; %to get rid of distweb(i,i)
Path=mean(distweb(distweb~=0 & distweb<=spe)); %infinite distances excluded in case of a not fully connected network

%% 12. TL: trophic level of each species
% a. T1: shortest path to a basal species
% T1=100*ones(spe,1); %infinite value
% tmp=false(spe,1);
% tmp(basal)=true;
% k=1;
% T1(tmp)=k;
% while any(T1==100) && k<=spe
%     k=k+1;
%     tmp=logical(sum(nicheweb(:,tmp)==1,2));
%     T1(tmp & T1==100)=k;
% end

%b. T2: Levine method. /!\ no cannibalism
W=zeros(spe);
nocanweb=nicheweb;
nocanweb(logical(eye(spe)))=0;
prey=sum(nocanweb,2)*ones(1,spe);
W(prey~=0)=nocanweb(prey~=0)./prey(prey~=0); %W_i_j: 1/number of preys of i if i eats j; 0 otherwise
T2=(inv(eye(spe)-W))*ones(spe,1);
if any(isinf(T2) | T2>spe+1)
    T2=repelem(NaN,spe)';
end
T=T2; %T=(T1+T2)/2
TL=mean(T);

%% 13. Omn: % of omnivorous species (eat at least 2 species of different trophic level)
TLprey=ones(spe,1)*round(T)'.*nicheweb; %row i contains the trophic level of each prey of i
omni=false(spe,1);
for i=1:spe
    TLpreyi=TLprey(i,:);
    omni(i)=(length(unique(TLpreyi(TLpreyi~=0)))>1); %omni(i)=true if the preys of i have at least 2 different trophic level
end
Omn=nnz(omni)/spe;

%% 14. Loop: % of species in loops (cannibalism excluded)
loop=false(spe,1); %loop(i)=true : i is part of at least one loop
tmp=double(nocanweb);
noloopweb=double(nocanweb);
for k=1:(nnz(nonbasal)-1) %we test path with k intermediary nodes (kmax : number of nonbasal species - 1)
    tmp=tmp*double(nocanweb); %tmp(i,j): number of path from predator i to prey j with k intermediary nodes
    loop_k=diag(tmp)>0; %tmp(i,i)>0 : there is at least one path from i to i with k indermediaty nodes
    loop(loop_k)=true;
    noloopweb(loop_k,loop_k)=0; %we eliminate all the links between loop_k species
end
Loop=nnz(loop)/spe;

%% 15. ChLen: average length of all food chains 
tmp=noloopweb;
chainnum=zeros(nnz(nonbasal),1);
for k=1:nnz(nonbasal) %/!\ k : number of links of the chain
    chainnum(k)=sum(sum(tmp(nonbasal,basal)),2);
    tmp=tmp*noloopweb;
end
chainlen=repelem(1:nnz(nonbasal),chainnum);
ChLen=mean(chainlen);

%% 16.ChNum: logarithm of the number of food chains 
ChNum=log(sum(chainnum));

%% 17. ChSD: standard deviation of all food chains
ChSD=std(chainlen);

%% 18. Con: connectance
Con=nnz(nicheweb)/spe^2;

%% 14bis, 15bis, 16bis, 17bis : Method2 : wrong method, can still imply loops.
% tmp=nocanweb; %tmp(i,j): number of path from predator i to prey j with k intermediary nodes
% for k=1:nnz(nonbasal) %/!\ k : number of links of the chain
%     chainnum(k)=sum(sum(tmp(nonbasal,basal)),2); %number of chain from predator i to basal species j with k links
%     loop_k=diag(tmp)>0; %tmp(i,i)>0 : there is at least one path from i to i with k links
%     loop(loop_k)=true;
%     tmp(logical(eye(spe)))=0; %elimination of the loops
%     tmp=tmp*nocanweb;
% end
% Loop=nnz(loop)/spe;
% chainlen=repelem(1:nnz(nonbasal),chainnum);
% ChLen=mean(chainlen);
% ChNum=log(sum(chainnum));
% ChSD=std(chainlen);

%% 
properties=[Bas, Top, Int, Can, Herb, MaxSim, VulSD, GenSD, ...
    LinkSD, Clust, Path, TL, Omn, Loop, ChLen, ChNum, ChSD, Con];

