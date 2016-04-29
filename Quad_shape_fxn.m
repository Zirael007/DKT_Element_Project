function N =  Quad_shape_fxn(xi, yi)

    % Shape functions for Quad element
    % 8 noded element (4 corners and 4 mid nodes)

    N(1) = -0.25*((1-xi)*(1-yi)*(1+xi+yi));
    N(2) = -0.25*((1+xi)*(1-yi)*(1-xi+yi));
    N(3) = -0.25*((1+xi)*(1+yi)*(1-xi-yi));
    N(4) = -0.25*((1-xi)*(1+yi)*(1+xi-yi));
    N(5) = 0.5*(1-xi^2)*(1-yi);
    N(6) = 0.5*(1+xi)*(1-yi^2);
    N(7) = 0.5*(1-xi^2)*(1+yi);
    N(8) = 0.5*(1-xi^2)*(1-yi^2);
end