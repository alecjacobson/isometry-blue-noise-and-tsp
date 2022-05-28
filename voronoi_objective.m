function [F,dFdX,I,D,C] = voronoi_objective(X,Y,rho)
  [I,D] = knnsearch(X,Y,'K',1);
  F = sum(rho.*D.^2);
  if nargout > 1
    n = size(X,1);
    % Voronoi masses 
    M = accumarray(I,rho,[n 1]);
    % Voronoi centers
    C = [accumarray(I,rho.*Y(:,1),[n 1]) accumarray(I,rho.*Y(:,2),[n 1])]./M;
    dFdX = 2*M.*(X-C);
  end
end
