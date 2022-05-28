X = readDMAT('secord-X-4096.dmat');
rng(0);
tic;
%[P,cost] = tsp(X,'spiral');
[P,cost] = tsp(X,'greedy-insertion');
toc
cost
plt(X(P,:),'-ok');
axis equal;
set(gca,'Ydir','reverse');
OP = P;
P = OP;

while true
I = edge_intersections(X,[P(1:end) P([2:end 1])]);
if isempty(I)
  break;
end
a = I(1,1);
c = I(1,2);
if a>c
  [c,a] = deal(a,c);
end
b = mod(a+1-1,numel(P))+1;
d = mod(c+1-1,numel(P))+1;

P0 = P;
if d==1
P = P([1:a c:-1:b ]);
else
P = P([1:a c:-1:b d:end]);
end
clf;
hold on;
%plt(X(P0([d:end 1:a]),:),'-ok');
plt(X(P,:),'-ok');

plt(X(P0([a b]),:),'-or','LineWidth',3);
plt(X(P0(b:c),:),'-og'  ,'LineWidth',3);
plt(X(P0([c d]),:),'-ob','LineWidth',3);
hold off;
axis equal;
set(gca,'Ydir','reverse');
drawnow;
assert(all(accumarray(P,1)==1))
        figgif('tsp.gif');

end
plt(X(P([1:end 1]),:),'-ok');
axis equal;
set(gca,'Ydir','reverse');


%tic;
%E = [P(1:end) P([2:end 1])];
%[B1,B2] = box_each_element(X,E);
%I = box_intersect(B1,B2);
%% prune incident edges
%keep = ...
%  E(I(:,1),1) ~= E(I(:,2),1) &  ...
%  E(I(:,1),1) ~= E(I(:,2),2) &  ...
%  E(I(:,1),2) ~= E(I(:,2),1) &  ...
%  E(I(:,1),2) ~= E(I(:,2),2);
%I = I(keep,:);
%sqrD = segment_segment_squared_distance( ...
%  X(E(I(:,1),1),:), X(E(I(:,1),2),:), X(E(I(:,2),1),:), X(E(I(:,2),2),:));
%I = I(sqrD<1e-7,:);
%toc


%    v--\       v---\
%   a    d      a    d
%    \  ^       |    ^
%     \/        |    |
%     /\        |    |
%    /  v       v    |
%   c    b      c    b
%    ^--/        \--^

%% This is bad idea. Lots of islands.
%% unsafe "symmetric" 2-opts
%C = normrow(diff(X(P,:)));
%A = sparse(P(1:end-1),P(2:end),C,size(X,1),size(X,1));
%%A = A+A';
%[I,J,V] = find(A);
%valid = true(size(I));
%BC = 0.5*(X(I,:)+X(J,:));
%[K,D] = knnsearch(BC,BC,'K',2);
%[~,S] = sort(D(:,2));
%for e1 = reshape(S,1,[])
%  e2 = K(e1,2);
%  if any(~valid([e1 e2]))
%    continue
%  end
%  if numel(unique([I(e1) J(e1) I(e2) J(e2)])) ~= 4
%    continue
%  end
%  clf;
%  hold on;
%  tsurf([I J],X);
%  plt(X([I(e1) J(e1)],:),'-or');
%  plt(X([I(e2) J(e2)],:),'-ob');
%  cur_cost = V(e1)+V(e2)
%  for pass = 1:2
%    switch pass
%    case 1
%      [i1,j1,i2,j2] = deal(I(e1),I(e2),J(e1),J(e2));
%    case 2
%      [i1,j1,i2,j2] = deal(I(e1),J(e2),J(e1),I(e2));
%    end
%    c1 = normrow(X(i1,:)-X(j1,:));
%    c2 = normrow(X(i2,:)-X(j2,:));
%    new_cost = c1 + c2;
%    if new_cost < cur_cost
%      valid([e1 e2]) = false;
%      [I(e1),J(e1),V(e1),I(e2),J(e2),V(e2)] = deal(i1,j1,c1,i2,j2,c2);
%      break;
%    end
%  end
%  
%  cur_cost
%  new_cost
%  plt(X([I(e1) J(e1)],:),'-oc');
%  plt(X([I(e2) J(e2)],:),'-og');
%  hold off;
%  axis equal;
%  set(gca,'Ydir','reverse');
%  drawnow
%end
