function B=DKQ_strain_displacement(J, dHxdxi, dHxdyi, dHydxi, dHydyi)

j = inv(J);

B=[j(1,1)*dHxdxi + j(1,2)*dHxdyi; 
    j(2,1)*dHydxi+j(2,2)*dHydyi;
    j(1,1)*dHydxi+j(1,2)*dHydyi+j(2,1)*dHxdxi+j(2,2)*dHxdyi];

end
