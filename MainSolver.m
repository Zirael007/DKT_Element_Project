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

load material.dat;
load geometry.dat;

l = geometry(:,1);
b = geometry(:,2);
t = geometry(:,3);

[coords, connectivity] = Mesh_Gen(l,b, [-0.5 -0.5], 20, 20);

nel = size(connectivity,1);
nnel = size(connectivity,2);
ndof = 3;
nnode = length(coords);
edof = ndof*nnel;
sdof = ndof*nnode;

PlotMesh(coords,connectivity);

force = zeros(sdof,1);
stif = zeros(sdof,sdof);
B = zeros(3,edof);


E = material(:,1);
nu = material(:,2);

P = 0.5;

D = (E*(t^3)/(12*(1-nu^2)))*[1 nu 0; nu 1 0; 0 0 0.5*(1-nu)];

%%

[pt,wt] = G_Quadrature('second');

% Loop over all elements

for iel = 1:1:nel
    
    % Extract coordinate data
    for i=1:nnel
        node(i)=connectivity(iel,i);
        x(i)=coords(node(i),1);
        y(i)=coords(node(i),2);
    end
    
    x_ij = -diff([x, x(1)]);
    y_ij = -diff([y, y(1)]);
    l_ij = x_ij.^2 + y_ij.^2;

    k = zeros(edof,edof); 
    f = zeros(edof,1);

    if nnel == 4
        [k, f] = DKQ_integration(x_ij, y_ij, l_ij, pt, wt);
    elseif nnel ==3
        [k, f] = DKT_integration(x_ij, y_ij, l_ij, pt, wt);
    end
    
    index = Elem_DOF(node,nnel,ndof);

    [stif, force] = Global_Assembly(stif, force, k, f, index);
    
end

bcdof = BoundaryCondition('ss-ss-ss-ss',coords) ;
bcval = zeros(1,length(bcdof)) ;

[stif,force] = constraints(stif,force,bcdof,bcval);

displacement = stif\force;

[w, thetax, thetay] = Post_Proc(displacement);