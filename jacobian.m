function J = jacobian(x_ij,y_ij, xi, yi)

J=0.25*[-x_ij(1)+x_ij(3)+yi*(x_ij(1)+x_ij(3)) -y_ij(1)+y_ij(3)+yi*(y_ij(1)+y_ij(3));
        -x_ij(2)+x_ij(4)+xi*(x_ij(1)+x_ij(3))  -y_ij(2)+y_ij(4)+xi*(y_ij(1)+y_ij(3))];

end