function N =  Quad_shape_fxn(x, y)

    % Shape functions for Quad element
    % 8 noded element (4 corners and 4 mid nodes)

    N(1) = *(1 - x -y)*(0.5 - x - y);
    N(2) = x*(2*x - 1);
    N(3) = y*(2*y - 1);
    N(4) = 4*x*y;
    N(5) = 4*y*(1 - x - y);
    N(6) = 4*x*(1 - x - y);
end