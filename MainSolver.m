%   Variable Descriptions:
%   nel     : No. of elements
%   nnel    : No. of nodes per element
%   ndof    : No. of DOF per node
%   nnode   : No. of nodes
%   edof    : No. of DOF per element
%   sdof    : No. of DOF for system
%   force   : System Force vector
%   stif    : System Stiffness matrix
%   B       : Element Strain-Displacement Matrix
%   l       : Length along X-axis
%   b       : Length along Y-axis
%   t       : Thickness of plate
%   E       : Young's Modulus
%   nu      : Poisson's Ratio
%
%%

clc
clear

load coords.dat;
load connectivity.dat;
load material.dat;
load geometry.dat;
nel = length(connectivity);
nnel = size(connectivity,2);
ndof = 3;
nnode = length(coords);
edof = ndof*nnel;
sdof = ndof*nnode;

force = zeros(sdof,1);
stif = zeros(sdof,sdof);
B = zeros(3,edof);

l = geometry(:,1);
b = geometry(:,2);
t = geometry(:,3);

E = material(:,1);
nu = material(:,2);

D = (E*t^3/(12*(1-nu)^2))*[1 nu 0; nu 1 0; 0 0 0.5*(1-nu)];

%%

% Loop over all elements

for iel = 1:1:nel
    
    % Extract coordinate data
    for i=1:nnel
        node(i)=connectivity(iel,i);
        x(i)=coords(node(i),1);
        y(i)=coords(node(i),2);
    end
    
end