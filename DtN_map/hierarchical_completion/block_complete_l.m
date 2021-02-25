function [block, mask, err] = block_complete_l(A, p, min_dim, r,tol)
% function to recursively compress an upper diagonal block

if min_dim >= size(A,1)
    
    block = A;
    mask = ones(size(A));
    err = 0.0;

else
    
    A11 = A(1:end/2,1:end/2);
    A12 = A(1:end/2,1+end/2:end);
    A21 = A(1+end/2:end,1:end/2);
    A22 = A(1+end/2:end,1+end/2:end);
    
    po = p;
    pd = min(1, p*r); 
    pl = min(1, p*r); 
    pu = min(1, p*r); 
    
    [block11, mask11, err11] = block_complete_o(A11,p ,min_dim, r,tol);
    [block12, mask12, err12] = block_complete_l(A12,pl,min_dim, r,tol);
    [block21, mask21, err21] = block_complete_o(A21,p ,min_dim, r,tol);
    [block22, mask22, err22] = block_complete_o(A22,p ,min_dim, r,tol);
    
    block = [block11, block12; ...
        block21, block22];
    
    mask = [mask11, mask12;...
        mask21, mask22];
    
    err = max([err11, err12, err21, err22]);
end


end