function mu = coherence(U)
    % functuon to compute the coherence of a subsaacep given by the columns
    % of the matrix U
    [n,r] = size(U);
    
    % we extract and orthonormal basis from U
    [Q, ~] = qr(U, 0);
    
    P = Q*(Q'*eye(n));
    
    normP = sum((P.^2),1);
    mu = n/r*max(normP);

end