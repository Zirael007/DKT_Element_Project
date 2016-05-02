    function [G_point,G_weight] = G_Quadrature(order)

switch order
    case 4 
   
        % 2-point quadrature rule
        G_point = [-0.577350269189626 -0.577350269189626;
                       0.577350269189626  0.577350269189626] ;
        G_weight = [1 1;1 1];        
    case 3
      
        % 3-point quadrature eule
        G_point = [0.5 0; 0.5 0.5; 0 0.5] ;
        G_weight = [1/3 1/3; 1/3 1/3; 1/3 1/3] ;
        
end

    end