
function B= DKT_strain_displacement(x_ij,y_ij,dHxdxi, dHxdyi, dHydxi, dHydyi)

A=x_ij(3)*y_ij(1)-x_ij(1)*y_ij(3);

B= 1/(2*A)*[ y_ij(3)*dHxdxi+y_ij(1)*dHxdyi;
            -x_ij(3)*dHydxi-x_ij(1)*dHydyi;
            -x_ij(3)*dHxdxi-x_ij(1)*dHxdyi+y_ij(3)*dHydxi+y_ij(1)*dHydyi];


end