function N =  Tri_shape_fxn(x, y)

    % Shape functions for Triangular element
    % 6 noded element (3 corners and 3 mid nodes)

    N(1) = 2*(1 - x -y)*(0.5 - x - y);
    N(2) = x*(2*x - 1);
    N(3) = y*(2*y - 1);
    N(4) = 4*x*y;
    N(5) = 4*y*(1 - x - y);
    N(6) = 4*x*(1 - x - y);

end