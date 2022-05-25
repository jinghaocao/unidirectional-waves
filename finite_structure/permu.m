function P = permu(x)
   n = length(x);
   M = eye(n);
   P = M;
   for j = 1:n
       P(j,:) = M(x(j),:);
   end
end