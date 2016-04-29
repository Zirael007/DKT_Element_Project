function B=DKQ_strain_displacement(J, dHxdxi, dHxdyi, dHydxi, dHydyi)


% j11 = J(2,2)/det(J);
% j12 = -J(1,2)/det(J);
% j21 = -J(2,1)/det(J);
% j22 = J(2,2)/det(J);

j = inv(J);

B=[j(1,1)*dHxdxi + j(1,2)*dHxdyi; 
    j(2,1)*dHydyi+j(2,2)*dHydxi;
    j(1,1)*dHydyi+j(1,2)*dHydxi+j(2,1)*dHxdxi+j(2,2)*dHxdyi];

end
