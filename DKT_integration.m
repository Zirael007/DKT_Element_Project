function [k, f] = DKT_integration(x_ij, y_ij, l_ij, pt, wt)

    % Loop over all integration points
       
       
    for intx = 1:1:2
        xi = pt(intx,1);
        wtx = wt(intx,1);
        for inty=1:1:2
            yi = pt(inty,2);
            wty = wt(inty,2);

            [dHxdxi,dHxdyi, dHydxi, dHydyi] = DKT_dH(p_k, q_k, r_k, t_k, xi, yi);
            
            N = DKT_shape_fxn(xi,yi);
            
            B = DKT_strain_displacement(x_ij, y_ij, dHxdxi, dHxdyi, dHydxi, dHydyi);
            
            k = k + B'*D*B*wtx*wty*det(J);
            
            fe = El_Force(nnel, N, P);
            f = f + fe*wtx*wty*det(J);
        end
    end

end