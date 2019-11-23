%% isConnected.m

%% Previous authors -------------------------------------------------------
% Perrine Tonin

%% Last update ------------------------------------------------------------
% when: 2-26-2018
% who: Valentin Cocco
% mail: valentin.cocco@ens.fr
% what : personal rewriting

%% Description ------------------------------------------------------------
% Check if one net work is fully connected (no independant networks)
% Inputs:
%   - M: matrix of adjacency ; M(i,j)>0 if i eats j. Otherwise M(i,j)=0
% Outputs:
%   - z: logical ; z=1 if fully connected network. 0 otherwise.
% Called by:
%   - NicheModel.m

%% 
function z=isConnected(M)
m=length(M);
M=max(M,M'); %undirected links: M(i,j)=1 if i eats j OR i is eaten by j (symetric matrix)
C=false(1,m); %C contains nodes connected amongst themselves (initiation : none)
N=[true, false(1,m-1)]; %N contains the new nodes connected to the nodes in C (initiation: node 1)
while any(N) %as long as we add new nodes
    C=C|N; %add the new nodes from N to C
    N=sum(M(N,:),1) & ~C; %select the nodes not already included in which are connected to the new nodes of C
end
z=all(C);
