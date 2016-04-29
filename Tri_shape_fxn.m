function N =  Tri_shape_fxn(xi, yi)

    % Shape functions for Triangular element
    % 6 noded element (3 corners and 3 mid nodes)

    N(1) = 2*(1 - xi -yi)*(0.5 - xi - yi);
    N(2) = xi*(2*xi - 1);
    N(3) = yi*(2*yi - 1);
    N(4) = 4*xi*yi;
    N(5) = 4*yi*(1 - xi - yi);
    N(6) = 4*xi*(1 - xi - yi);

end