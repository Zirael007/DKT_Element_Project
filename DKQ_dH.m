function [dHxdxi,dHxdyi, dHydxi, dHydyi]=DKQ_dH(coords,xi,yi)
  
    
%            x12=coords(1,1)-coords(2,1);
%            x23=coords(2,1)-coords(3,1);
%            x34=coords(3,1)-coords(4,1);
%            x41=coords(4,1)-coords(1,1);
%            y12=coords(1,2)-coords(2,2);
%            y23=coords(2,2)-coords(3,2);
%            y34=coords(3,2)-coords(4,2);
%            y41=coords(4,2)-coords(1,2);
           rot_coords = [coords; coords(1,:)];

           x_ij = diff(rot_coords(:,1));
           y_ij = diff(rot_coords(:,2));
           
           l_ij = x_ij.^2 + y_ij.^2;
           
%            l_sq12=x12^2+y12^2;
%            l_sq23=x23^2+y23^2;
%            l_sq34=x34^2+y34^2;
%            l_sq41=x41^2+y41^2;

           a_k = num2cell(-x_ij./l_ij);
           [a5, a6, a7, a8] = a_k{:};
           

%            a5=-x12/l_sq12;
%            a6=-x23/l_sq23;
%            a7=-x34/l_sq34;
%            a8=-x41/l_sq41;

           b_k = num2cell(0.75.*x_ij.*y_ij./l_ij);
           [b5, b6, b7, b8] = b_k{:};

%            b5=0.75*x12*y12/l_sq12;
%            b6=0.75*x23*y23/l_sq12;
%            b7=0.75*x12*y34/l_sq34;
%            b8=0.75*x12*y41/l_sq41;

           c_k = num2cell((0.25.*x_ij.^2 - 0.5.*y_ij.^2)./l_ij);
           [c5, c6, c7, c8] = c_k{:};

%            c5=(0.25*x12^2-0.5*y12^2)/l_sq12;
%            c6=(0.25*x23^2-0.5*y23^2)/l_sq23;
%            c7=(0.25*x34^2-0.5*y34^2)/l_sq34;
%            c8=(0.25*x41^2-0.5*y41^2)/l_sq41;

           d_k = num2cell(-y_ij.^2./l_ij);
           [d5, d6, d7, d8] = d_k{:};

%            d5=-y12^2/l_sq12;
%            d6=-y23^2/l_sq23;
%            d7=-y34^2/l_sq34;
%            d8=-y41^2/l_sq41;

           e_k = num2cell((0.25.*y_ij.^2 - 0.5.*x_ij.^2)./l_ij);
           [e5, e6, e7, e8] = e_k{:};

%            e5=(0.25*y12^2-0.5x12^2)/l_sq12;
%            e6=(0.25*y23^2-0.5x23^2)/l_sq23;
%            e7=(0.25*y34^2-0.5x34^2)/l_sq34;
%            e8=(0.25*y41^2-0.5x41^2)/l_sq41;
           
           dNdxi=[0.25*(2*xi+yi)*(1-yi),...
                  0.25*(2*xi-yi)*(1-yi),...
                  0.25*(2*xi+yi)*(1+yi),...
                  0.25*(2*xi-yi)*(1+yi),...
                  -xi*(1-yi),...
                  0.5*(1-yi^2),...
                  -xi*(1+yi),...
                  -0.5*(1-yi^2)];
           
           dNdyi=[0.25*(2*yi+xi)*(1-yi),...
                  0.25*(2*yi-xi)*(1+yi),...
                  0.25*(2*yi+xi)*(1+yi),...
                  0.25*(2*yi-xi)*(1-yi),...
                  -0.5*(1-yi^2),...
                  -yi*(1+xi),...
                  0.5*(1-yi^2),...
                  -yi*(1-xi)];
          
          dHxdxi=[1.5*(a5*dNdxi(5)-a8*dNdxi(8));
                   b5*dNdxi(5)+b8*dNdxi(8);
                   dNdxi(1)-c5*dNdxi(5)-c8*dNdxi(8);
                   1.5*(a6*dNdxi(6)-a5*dNdxi(5));
                   b6*dNdxi(6)+b5*dNdxi(5);
                   dNdxi(2)-c6*dNdxi(6)-c5*dNdxi(5);
                   1.5*(a7*dNdxi(7)-a6*dNdxi(6));
                   b7*dNdxi(7)+b6*dNdxi(6);
                   dNdxi(3)-c7*dNdxi(7)-c6*dNdxi(6);
                   1.5*(a8*dNdxi(8)-a7*dNdxi(7));
                   b8*dNdxi(8)+b7*dNdxi(7);
                   dNdxi(4)-c8*dNdxi(8)-c7*dNdxi(7)];
                  
                   
          dHxdyi=[1.5*(a5*dNdyi(5)-a8*dNdyi(8));
                   b5*dNdyi(5)+b8*dNdyi(8);
                   dNdyi(1)-c5*dNdyi(5)-c8*dNdyi(8);
                   1.5*(a6*dNdyi(6)-a5*dNdyi(5));
                   b6*dNdyi(6)+b5*dNdyi(5);
                   dNdyi(2)-c6*dNdyi(6)-c5*dNdyi(5);
                   1.5*(a7*dNdyi(7)-a6*dNdyi(6));
                   b7*dNdyi(7)+b6*dNdyi(6);
                   dNdyi(3)-c7*dNdyi(7)-c6*dNdyi(6);
                   1.5*(a8*dNdyi(8)-a7*dNdyi(7));
                   b8*dNdyi(8)+b7*dNdyi(7);
                   dNdyi(4)-c8*dNdyi(8)-c7*dNdyi(7)];
                   
          
          dHydxi=[1.5*(d5*dNdxi(5)-d8*dNdxi(8));
                   -dNdxi(1)-e5*dNdxi(5)-e8*dNdxi(8);
                   -b5*dNdxi(5)-b8*dNdxi(8);
                   1.5*(d7*dNdxi(7)-d6*dNdxi(6));
                   -dNdxi(2)-e6*dNdxi(6)-e5*dNdxi(5);
                   -b6*dNdxi(6)-b5*dNdxi(5);
                   1.5*(d7*dNdxi(7)-d6*dNdxi(6));
                   -dNdxi(3)-e7*dNdxi(7)-e6*dNdxi(6);
                   -b7*dNdxi(7)-b6*dNdxi(6);
                   1.5*(d8*dNdxi(8)-d7*dNdxi(7))
                   -dNdxi(4)-e8*dNdxi(8)-e7*dNdxi(7);
                   -b8*dNdxi(8)-b7*dNdxi(7)];
               
           dHydyi=[1.5*(d5*dNdyi(5)-d8*dNdyi(8));
                   -dNdyi(1)-e5*dNdyi(5)-e8*dNdyi(8);
                   -b5*dNdyi(5)-b8*dNdyi(8);
                   1.5*(d7*dNdyi(7)-d6*dNdyi(6));
                   -dNdyi(2)-e6*dNdyi(6)-e5*dNdyi(5);
                   -b6*dNdyi(6)-b5*dNdyi(5);
                   1.5*(d7*dNdyi(7)-d6*dNdyi(6));
                   -dNdyi(3)-e7*dNdyi(7)-e6*dNdyi(6);
                   -b7*dNdyi(7)-b6*dNdyi(6);
                   1.5*(d8*dNdyi(8)-d7*dNdyi(7))
                   -dNdyi(4)-e8*dNdyi(8)-e7*dNdyi(7);
                   -b8*dNdyi(8)-b7*dNdyi(7)];
end