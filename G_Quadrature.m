    function [G_point,G_weight] = G_Quadrature(order)
%-------------------------------------------------------------------
% Depending on the order, this function determines the Gauss Quadrature
% based on 2- point or 1-point quadrature rules. 
%-------------------------------------------------------------------

switch order
    case 'second' 
   
        % 2-point quadrature rule
        G_point = [-0.577350269189626 -0.577350269189626;
                       0.577350269189626  0.577350269189626] ;
        G_weight = [1 1;1 1];        
    case 'first'
      
        % 1-point quadrature eule
        G_point = [0 0] ;
        G_weight = [2. 2.] ;
        
end

    end