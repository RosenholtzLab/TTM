function out = getBiggestNonNanSquare(A, maxSizeOfA)
    [ht wt] = size(A);
    assert(ht == wt);
    mid = (ht + 1)/2;
    for span = (mid-1) : -1 : 0
        out = A((-span:span)+mid, (-span:span)+mid);
        if (any(isnan(out(:))) | size(out,1)>maxSizeOfA)
            continue;
        else
            return;
        end
    end
    out = [];
return