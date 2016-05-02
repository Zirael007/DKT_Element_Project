function N =  DKQ_shape_fxn(xi, yi)

    % Shape functions for Quad element
    % 4 noded element (Kirchoff Plate Element)

    N(1) = -0.25*((1-xi)*(1-yi)*(1+xi+yi));
    N(2) = -0.25*((1+xi)*(1-yi)*(1-xi+yi));
    N(3) = -0.25*((1+xi)*(1+yi)*(1-xi-yi));
    N(4) = -0.25*((1-xi)*(1+yi)*(1+xi-yi));
    N(5) = 0.5*(1-xi^2)*(1-yi);
    N(6) = 0.5*(1+xi)*(1-yi^2);
    N(7) = 0.5*(1-xi^2)*(1+yi);
    N(8) = 0.5*(1-xi^2)*(1-yi^2);
end