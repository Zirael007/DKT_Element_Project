function PlotMesh(coords,connectivity)

nel = size(connectivity,1) ;                  % number of elements
nnode = size(coords,1) ;          % total number of connectivity in system
nnel = size(connectivity,2);                % number of connectivity per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;

for iel=1:nel   
     for i=1:nnel
     nd(i)=connectivity(iel,i);         % extract connected node for (iel)-th element
     X(i,iel)=coords(nd(i),1);    % extract x value of the node
     Y(i,iel)=coords(nd(i),2);    % extract y value of the node
     end
end
    
% Plotting the FEM mesh, diaplay Node numbers and Element numbers
     f1 = figure ;
     set(f1,'name','Mesh','numbertitle','off','Color','w') ;
     fill(X,Y,'w')
     
     title('Finite Element Mesh') ;
     axis off ;
     
% To disply the node numbers     
     k = connectivity(:,1:end);
     nd = k' ;
    for i = 1:nel
        text(X(:,i),Y(:,i),int2str(nd(:,i)),'fontsize',8,'color','k');
        text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
    end
end