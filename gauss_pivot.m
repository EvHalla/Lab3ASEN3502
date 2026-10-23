function x = gauss_pivot(A, b)
n = size(b,1);
aug = [A,b];
for c = 1:n
    [~, best_idx] = max(abs(aug([c:end],c)));
    best_idx = best_idx+c-1;
    aug([best_idx, c],:) = aug([c,best_idx],:); 
    for r = c+1:n
        f = aug(r,c)/aug(c,c);
        aug(r, c:end) = aug(r, c:end) - f * aug(c,c:end);
    end
end

x = zeros(n,1);
x(n) = aug(n,end) / aug(n,n);
for i = n-1:-1:1
    x(i) = (aug(i, end) - aug(i,i+1:n) * x(i+1:n)) / aug(i,i);
end

return


