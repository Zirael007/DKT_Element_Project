function [kk,ff]=Global_Assembly(kk,ff,k,f,index)
 
 ele_dof = length(index);
 for i=1:ele_dof
   ii=index(i);
     ff(ii)=ff(ii)+f(i);
     for j=1:ele_dof
       jj=index(j);
         kk(ii,jj)=kk(ii,jj)+k(i,j);
     end
 end

end