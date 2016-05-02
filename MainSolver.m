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

[coords, connectivity] = Mesh_Gen(l,b, [-0.5 -0.5], 10, 10, 'quad');

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

P = 5000;

D = (E*(t^3)/(12*(1-nu^2)))*[1 nu 0; nu 1 0; 0 0 0.5*(1-nu)];

%%

[pt,wt] = G_Quadrature(nnel);

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
    
    
    a_k = num2cell(-x_ij./l_ij);
    b_k = num2cell(0.75.*x_ij.*y_ij./l_ij);
    c_k = num2cell((0.25.*x_ij.^2 - 0.5.*y_ij.^2)./l_ij);
    d_k = num2cell(-y_ij.^2./l_ij);
    e_k = num2cell((0.25.*y_ij.^2 - 0.5.*x_ij.^2)./l_ij);

    p_k = num2cell(-6*x_ij./l_ij);
    q_k = num2cell(3*x_ij.*y_ij./l_ij);
    r_k = num2cell(3*y_ij.^2./l_ij);
    t_k = num2cell(-6*y_ij./l_ij);       
    
    k = zeros(edof,edof); 
    f = zeros(edof,1);

    % Loop over all integration points
        
    k = zeros(edof,edof); 
    f = zeros(edof,1);
    

    
    if nnel == 4
%        [k, f] = DKQ_integration(x_ij, y_ij, l_ij, pt, wt, nnel, edof, D, P);
        for intx = 1:1:2
        xi = pt(intx,1);
        wtx = wt(intx,1);
            for inty=1:1:2
                yi = pt(inty,2);
                wty = wt(inty,2);

                [dHxdxi,dHxdyi, dHydxi, dHydyi] = DKQ_dH(a_k, b_k, c_k, d_k, e_k, xi, yi);
            
                J = DKQ_jacob(x_ij, y_ij, xi, yi);
            
                N = DKQ_shape_fxn(xi,yi);
            
                B = DKQ_strain_displacement(J, dHxdxi, dHxdyi, dHydxi, dHydyi);
            
                k = k + B'*D*B*wtx*wty*det(J);
            
                fe = El_Force(nnel, N, P);
                f = f + fe*wtx*wty*det(J);
            end
        end
    elseif nnel ==3
%         [k, f] = DKT_integration(x_ij, y_ij, l_ij, pt, wt);
        for intx = 1:1:3
        xi = pt(intx,1);
        wtx = wt(intx,1);
         for inty=1:1:3
                yi = pt(inty,2);
                wty = wt(inty,2);

                [dHxdxi, dHxdyi, dHydxi, dHydyi] = DKT_dH(p_k, q_k, r_k, t_k, xi, yi);
            
%                 J = DKQ_jacob(x_ij, y_ij, xi, yi);
                A = x_ij(3)*y_ij(1)-x_ij(1)*y_ij(3);
            
                N = DKT_shape_fxn(xi,yi);
            
                B = DKT_strain_displacement(x_ij,y_ij,dHxdxi, dHxdyi, dHydxi, dHydyi);
            
                k = k + A*B'*D*B*wtx*wty;
            
                fe = El_Force(nnel, N, P);
                f = f + fe*wtx*wty*A;
          end
        end
        
    end
    
    index = Elem_DOF(node,nnel,ndof);

    [stif, force] = Global_Assembly(stif, force, k, f, index);
    
end

bcdof = BoundaryCondition('ss-ss-ss-ss',coords) ;
bcval = zeros(1,length(bcdof)) ;

[stif,force] = constraints(stif,force,bcdof,bcval);

displacement = stif\force;

[w, thetax, thetay] = Post_Proc(displacement);

fprintf('Max Tranverse Displacement: %f\n',min(w));