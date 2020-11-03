function   [Fmatrix BasisNum] = BasisSearch(AoA_range, Nt, epsilon)


A = [];
for ix = 1:1:Nt
    for iy = 1:1:Nt
        if ix == iy
            A(ix, iy) = max(AoA_range) - min(AoA_range);
        else
            A(ix, iy) = 1/(j*(ix-iy)*pi)*(exp(j*(ix-iy)*max(AoA_range)*pi) - exp(j*(ix-iy)*min(AoA_range)*pi));
        end
    end
end 

[ss vv dd] = svd(A);

Sum_Power = sum(diag(vv*vv));

ratio = diag(vv*vv)/Sum_Power;

BasisNum = 1;
while sum(ratio(1:BasisNum))<epsilon
    BasisNum = BasisNum + 1;
end

Fmatrix = dd(:,1:BasisNum);