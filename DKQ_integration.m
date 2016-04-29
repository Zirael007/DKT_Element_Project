function [k, f] = DKQ_integration(x_ij, y_ij, l_ij, pt, wt, edof, nnel, D, P)

    % Loop over all integration points
    
    a_k = num2cell(-x_ij./l_ij);
    b_k = num2cell(0.75.*x_ij.*y_ij./l_ij);
    c_k = num2cell((0.25.*x_ij.^2 - 0.5.*y_ij.^2)./l_ij);
    d_k = num2cell(-y_ij.^2./l_ij);
    e_k = num2cell((0.25.*y_ij.^2 - 0.5.*x_ij.^2)./l_ij);

    k = zeros(edof,edof); 
    f = zeros(edof,1);
    
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

end