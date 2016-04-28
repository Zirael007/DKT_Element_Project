function [H_x,H_y]=H_derivatives(coords,eps,eta)
  
    
           x12=coords(1,1)-coords(2,1);
           x23=coords(2,1)-coords(3,1);
           x34=coords(3,1)-coords(4,1);
           x41=coords(4,1)-coords(1,1);
           y12=coords(1,2)-coords(2,2);
           y23=coords(2,2)-coords(3,2);
           y34=coords(3,2)-coords(4,2);
           y41=coords(4,2)-coords(1,2);
           
           l_sq12=x12^2+y12^2;
           l_sq23=x23^2+y23^2;
           l_sq34=x34^2+y34^2;
           l_sq41=x41^2+y41^2;
           
           a5=-x12/l_sq12;
           a6=-x23/l_sq23;
           a7=-x34/l_sq34;
           a8=-x41/l_sq41;
           
           b5=0.75*x12*y12/l_sq12;
           b6=0.75*x23*y23/l_sq12;
           b7=0.75*x12*y34/l_sq34;
           b8=0.75*x12*y41/l_sq41;
           
           c5=(0.25*x12^2-0.5*y12^2)/l_sq12;
           c6=(0.25*x23^2-0.5*y23^2)/l_sq23;
           c7=(0.25*x34^2-0.5*y34^2)/l_sq34;
           c8=(0.25*x41^2-0.5*y41^2)/l_sq41;
           
           d5=-y12^2/l_sq12;
           d6=-y23^2/l_sq23;
           d7=-y34^2/l_sq34;
           d8=-y41^2/l_sq41;
           
           e5=(0.25*y12^2-0.5x12^2)/l_sq12;
           e6=(0.25*y23^2-0.5x23^2)/l_sq23;
           e7=(0.25*y34^2-0.5x34^2)/l_sq34;
           e8=(0.25*y41^2-0.5x41^2)/l_sq41;   
           
           N_eps=[0.25*(2*eps+eta)*(1-eta);
                  0.25*(2*eps-eta)*(1-eta);
                  0.25*(2*eps+eta)*(1+eta);
                  0.25*(2*eps-eta)*(1+eta);
                  -eps*(1-eta);
                  0.5*(1-eta^2);
                  -eps*(1+eta);
                  -0.5*(1-eta^2)];
           
           N_eta=[0.25*(2*eta+eps)*(1-eta);
                  0.25*(2*eta-eps)*(1+eta);
                  0.25*(2*eta+eps)*(1+eta);
                  0.25*(2*eta-eps)*(1-eta);
                  -0.5*(1-eta^2);
                  -eta*(1+eps);
                  0.5*(1-eta^2);
                  -eta*(1-eps)];
          
          H_x_eps=[1.5*(a5*N_eps(5,1)-a8*N_eps(8,1);
                   b5*N_eps(5,1)+b8*N_eps(8,1);
                   N_eps(1,1)-c5*N_eps(5,1)-c8*N_eps(8,1);
                   1.5*(a6*N_eps(6,1)-a5*N_eps(5,1));
                   b6*N_eps(6,1)+b5*N_eps(5,1);
                   N_eps(2,1)-c6*N_eps(6,1)-c5*N_eps(5,1);
                   1.5*(a7*N_eps(7,1)-a6*N_eps(6,1));
                   b7*N_eps(7,1)+b6*N_eps(6,1);
                   N_eps(3,1)-c7*N_eps(7,1)-c6*N_eps(6,1);
                   1.5*(a8*N_eps(8,1)-a7*N_eps(7,1));
                   b8*N_eps(8,1)+b7*N_eps(7,1);
                   N_eps(4,1)-c8*N_eps(8,1)-c7*N_eps(7,1)];
                  
                   
          H_x_eta=[1.5*(a5*N_eta(5,1)-a8*N_eta(8,1);
                   b5*N_eta(5,1)+b8*N_eta(8,1);
                   N_eta(1,1)-c5*N_eta(5,1)-c8*N_eta(8,1);
                   1.5*(a6*N_eta(6,1)-a5*N_eta(5,1));
                   b6*N_eta(6,1)+b5*N_eta(5,1);
                   N_eta(2,1)-c6*N_eta(6,1)-c5*N_eta(5,1);
                   1.5*(a7*N_eta(7,1)-a6*N_eta(6,1));
                   b7*N_eta(7,1)+b6*N_eta(6,1);
                   N_eta(3,1)-c7*N_eta(7,1)-c6*N_eta(6,1);
                   1.5*(a8*N_eta(8,1)-a7*N_eta(7,1));
                   b8*N_eta(8,1)+b7*N_eta(7,1);
                   N_eta(4,1)-c8*N_eta(8,1)-c7*N_eta(7,1)];
                   
          
          H_y_eps=[1.5*(d5*N_eps(5,1)-d8*N_eps(8,1));
                   -N_eps(1,1)-e5*N_eps(5,1)-e8*N_eps(8,1);
                   -b5*N_eps(5,1)-b8*N_eps(8,1);
                   1.5*(d7*N_eps(7,1)-d6*N_eps(6,1));
                   -N_eps(2,1)-e6*N_eps(6,1)-e5*N_eps(5,1);
                   -b6*N_eps(6,1)-b5*N_eps(5,1);
                   1.5*(d7*N_eps(7,1)-d6*N_eps(6,1));
                   -N_eps(3,1)-e7*N_eps(7,1)-e6*N_eps(6,1);
                   -b7*N_eps(7,1)-b6*N_eps(6,1);
                   1.5*(d8*N_eps(8,1)-d7*N_eps(7,1))
                   -N_eps(4,1)-e8*N_eps(8,1)-e7*N_eps(7,1);
                   -b8*N_eps(8,1)-b7*N_eps(7,1)];
               
           H_y_eta=[1.5*(d5*N_eta(5,1)-d8*N_eta(8,1));
                   -N_eta(1,1)-e5*N_eta(5,1)-e8*N_eta(8,1);
                   -b5*N_eta(5,1)-b8*N_eta(8,1);
                   1.5*(d7*N_eta(7,1)-d6*N_eta(6,1));
                   -N_eta(2,1)-e6*N_eta(6,1)-e5*N_eta(5,1);
                   -b6*N_eta(6,1)-b5*N_eta(5,1);
                   1.5*(d7*N_eta(7,1)-d6*N_eta(6,1));
                   -N_eta(3,1)-e7*N_eta(7,1)-e6*N_eta(6,1);
                   -b7*N_eta(7,1)-b6*N_eta(6,1);
                   1.5*(d8*N_eta(8,1)-d7*N_eta(7,1))
                   -N_eta(4,1)-e8*N_eta(8,1)-e7*N_eta(7,1);
                   -b8*N_eta(8,1)-b7*N_eta(7,1)];
            
            H_x(1,:)=H_x_eps;
            H_x(2,:)=H_x_eta;
            H_y(1,:)=H_y_eps;
            H_y(2,:)=H_y_eta;
            
           
         
          
           
   




end