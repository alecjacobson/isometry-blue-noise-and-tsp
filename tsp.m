function [P,cost] = tsp(X,method)




  costfun = @(P) sum(normrow(diff(X(P([1:end 1]),:))));
  cost = [];
  P = [];

  switch method
  case 'spiral'
    [V,F] = spiral_mesh(ceil(sqrt(size(X,1))/2));
    V = (V-min(V))./(max(V)-min(V)).*(max(X)-min(X))+min(X);
    [B,L] = ordered_outline(F);
    len = normrow(diff(V(B([1:end 1]),:)));
    T = cumsum([0;len]);
    E = [B(1:end)' B([2:end 1])'];
    assert(size(E,2) == 2);
    [~,J,C] = point_mesh_squared_distance(X,V,E);
    L = barycentric_coordinates(C,V(E(J,1),:),V(E(J,2),:));
    G = T(J) + L(:,2).*len(J);
    [~,P] = sort(G);
    %tsurf(E,V);
    %hold on;
    %sct(X,'ok');
    %hold off;
    %axis equal;

  case 'greedy-insertion'
    % Initial curve
    H = convhull(X);
    H = H(1:end-1);
    I = setdiff((1:size(X,1))',H);

    E = [H(1:end) H([2:end 1])]';
    X2E = sparse(E,repmat(1:size(E,2),2,1),1,size(X,1),size(X,1));
    while ~isempty(I)
      [sqrD,J] = point_mesh_squared_distance(X(I,:),X,E');

      [~,s] = min(sqrD);
      j = J(s);
      i = I(s);
      % swap with end and remove
      [I(end),I(s)] = deal(I(s),I(end));
      I(end) = [];
      on_vertex = 0;
      for c = 1:2
        v = E(c,j);
        sqrDv = sum((X(v,:) - X(i,:)).^2);
        if sqrDv < sqrD(s)+1e-10
          on_vertex = c;
          % Is j really best edge to use?
          aj = det([X(E(1,j),:)-X(i,:);X(E(2,j),:)-X(i,:)]);
          o = setdiff(find(X2E(v,:)),j);
          ao = det([X(E(1,o),:)-X(i,:);X(E(2,o),:)-X(i,:)]);
          if ao>0
            j = o;
          end
          break;
        end
      end

      s = E(1,j);
      d = E(2,j);
      X2E(E(2,j),j) = 0;
      E(2,j) = i;
      X2E(E(2,j),j) = 1;
      n = size(E,2)+1;
      E(:,n) = [i;d];
      X2E([i;d],n) = 1;

      if mod(numel(I),100) == 1
        clf;
        hold on;
        tsurf(E',X);
        sct(X(E(:),:),'ok');
        sct(X(I,:),'.r');
        sct(X(i,:),'ob','LineWidth',3);
        hold off;
        axis equal;
        set(gca,'Ydir','reverse');
        figgif('tsp.gif');
        drawnow
      end
      %if ~isempty(intersections(X,E'))
      %  clf;
      %  hold on;
      %  tsurf(E',X);
      %  sct(X(E(:),:),'ok');
      %  sct(X(I,:),'.r');
      %  sct(X(i,:),'ob','LineWidth',3);
      %  hold off;
      %  axis equal;
      %  set(gca,'Ydir','reverse');
      %  keyboard
      %end
    end
  case 'greedy-selection'
    D = pdist2(X,X);
    D(sub2ind(size(D),1:size(D,1),1:size(D,1))) = inf;
    [I,J,V] = find(triu(D,1));
    [~,S] = sort(V);
    I = I(S);
    J = J(S);
    E = zeros(2,0);
    % counts
    C = zeros(size(X,1),1);
    % component
    K = zeros(size(X,1),1);
    CC = sparse(size(X,1),size(X,1));

    k = 1;
    num_fully_visited = 0;
    for iter = 1:numel(I)
      i = I(iter);
      j = J(iter);
      % Would adding this edge create a non-manifold vertex?
      if C(i)==2 || C(j)==2
        continue;
      end

      % Adding a fresh edge:
      % #starts → #starts+1
      % #ends   → #ends  +1
      % Adding to existing chain:
      % #starts → #starts
      % #ends   → #ends  
      % merging two existing chains:
      % #starts → #starts-1
      % #ends   → #ends-1
      % existing chain to loop
      % #starts → #starts-1
      % #ends   → #ends-1


      % Would adding this edge create a loop?
      %new_num_fully_visited = num_fully_visited + C(i,2) + C(j,1);
      %if ...
      %  size(E,2) < size(X,1)-1  &&  ... loop ok at very end
      %  new_num_fully_visited == size(E,1)
      %  continue;
      %end

      %if iter == 654
      %  clf;
      %  hold on;
      %  tsurf(E',X);
      %  plt(X([i,j],:),'-or');
      %  hold off;
      %  axis equal;
      %  set(gca,'Ydir','reverse');
      %  title(sprintf('iter: %d',iter),'Fontsize',20);
      %  %axis([0.0057    0.3312    0.7896    0.9963]);
      %  %pause
      %end

      %
      % o-1-o
      % o-1-o                o-2-o
      % o-1-o-1-o            o-2-o
      

      %num_fully_visited = new_num_fully_visited;
      if K(i) == 0 && K(j) == 0
        % new component
        K(i) = k;
        K(j) = k;
        CC([i,j],k) = 1;
        k = k+1;
      elseif K(i) == 0 
        K(i) = K(j);
        CC(i,K(i)) = 1;
      elseif K(j) == 0 
        K(j) = K(i);
        CC(j,K(j)) = 1;
      elseif K(i) == K(j)
        continue;
      elseif K(i) ~=0 && K(j) ~=0
        %if iter == 654; keyboard; end
        if K(i) < K(j)
          mink = K(i);
          maxk = K(j);
        else
          mink = K(j);
          maxk = K(i);
        end
        CC(:,mink) = CC(:,mink) + CC(:,maxk) ;
        CC(:,maxk) = 0;
        K(find( CC(:,mink) )) = mink;
      end
      C(i) = C(i)+1;
      C(j) = C(j)+1;
      assert(K(i) == K(j));
      E(:,end+1) = [i;j];
      %cc = @(A) conncomp(A);
      %cc = @(A) cc(A(any(A,1),any(A,1)));
      %if ( numel(unique(K(K~=0))) ~= cc(adjacency_matrix(E')))
      %  clf;
      %  hold on;
      %  %tsurf(E',X);
      %  tsurf(E',X,'VertexIndices',1);
      %  plt(X([i,j],:),'-or');
      %  hold off;
      %  axis equal;
      %  set(gca,'Ydir','reverse');
      %  title(sprintf('iter: %d',iter),'Fontsize',20);
      %  axis([0.6356    0.7044    0.8849    1.0298]);
      %  %axis([0.0057    0.3312    0.7896    0.9963]);
      %  %pause
      %  keyboard
      %end
    end

  case 'nearest-neighbor'
    D = pdist2(X,X);
    D(sub2ind(size(D),1:size(D,1),1:size(D,1))) = inf;
    P = [1];
    M = ones(size(X,1),1);
    while numel(P)<size(X,1)-1
      last = P(end);
      M(last) = inf;
      [~,next] = min(M.*D(:,last));
      P(end+1,:) = next;
    end
  case 'snake'
    h = normrow(max(X)-min(X))/sqrt(size(X,1));
    [~,B] = histc(X(:,1),min(X(:,1))-h/2:h:max(X(:,1))+h/2);
    P = [];
    for b = 1:max(B)
      % can factor out this O(n) operation
      Ib = find(B==b);
      % can factor out this O(n log n) by sorting everybody by X(:,2)
      [~,Yb] = sort(X(Ib,2));
      if mod(b,2) == 1
        Yb = flip(Yb);
      end
      P = [P;Ib(Yb)];
    end
  case 'angular'
    X = X-mean(X);
    A = atan2(X(:,2),X(:,1));
    [~,P] = sort(A);
  case 'random'
    max_iter = ceil(sqrt(size(X,1)));
    P = [];
    cost = inf;
    for iter = 1:max_iter
      Pi = randperm(size(X,1));
      costi = costfun(Pi);
      if costi < cost
        P = Pi;
        cost = costi;
      end
    end
  end
  if isempty(P)
    sd = find(accumarray(E(:),1)==1);
    [P,L] = ordered_outline([E sd]');
    P = reshape(P,[],1);
    assert(numel(L)==2);
  end
  if isempty(cost)
    cost = costfun(P);
  end
end
