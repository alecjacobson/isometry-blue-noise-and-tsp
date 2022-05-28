function Y = add_swirls(X,r,stride)
  if isempty(r)
    [~,D] = knnsearch(X,X,'K',2);
    r = mean(D(:,2))/10;
  end

  T = normalizerow(X-X([end 1:end-1],:));
  N = [T(:,2) -T(:,1)];
  A = atan2(N(:,2),N(:,1));
  %A = A([2:end 1])-A;
  ease = @(t) 3.*t.^2-2.*t.^3;ease = @(t) ease(min(max(t,0),1));
  Y = zeros(2,0);
  for i = 1:size(X,1)
    n = mod(i+1-1,size(X,1))+1;
    p = mod(i-1-1,size(X,1))+1;
    % rotate 90° CW
    %Y(:,end+1) = (X(i,:)+r*N(i,:))';
    %Y(:,end+1) = (X(i,:)+r*[cos(A(n)) sin(A(n))])';
    Ai = A(i);
    An = A(n);
    if An<Ai
      An = An+2*pi;
    end
    l = 2;
    An = An+l*2*pi;
    assert(An>Ai);
    ri = (rand(1,l)*2-1)*0.4;
    for t = linspace(0,1,ceil((An-Ai)/0.2)+1)
      a = Ai+t*(An-Ai);
      rt = r*(1+interp1(linspace(0,1,l),ri,t)*sin(ease(t)*l*2*pi));
      Y(:,end+1) = (X(i,:)+rt*[cos(a) sin(a)])';
    end
  end
  Y = Y';

end
