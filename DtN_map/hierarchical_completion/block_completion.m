function [block, E, err] = block_completion(A, p, tol)
% function to compute the matrix completion from a bloack A, by sampling
% the matrix with a Bernouilli mask of parameter p
% output block: compressed block
%        E : the mask with the sampled entries
%        err: error in the reconstruction

% if the sampling probability is one, sample the 
if p == 1
    block = A;
    E = ones(size(A));
    err = 0;
else
    
    % size of the
    [n,m] = size(A);
    
    % computing the mask
    E = (rand(n,m) < p);
    
    % indices of the mask
    indE = find(E);
    fprintf(strcat("number of data points is ", num2str(length(indE))))
    
    % we start the optimiziation
    cvx_precision best
    cvx_begin
    variable Xe(n,m)
    minimize norm_nuc(Xe)
    subject to
    Xe(indE) == A(indE)
    cvx_end
    
    % we consider the \ell^{\infty} norm to make things easier across
    % different sizes
    err = max(abs((Xe(:)-A(:))));
    
    fprintf(['Error on the constraint :', num2str(norm(Xe(indE)-A(indE))), '\n'])
    fprintf(['Error on the reconstruction :', num2str(norm(Xe(:)-A(:))), '\n'])
    
    block = Xe;
    
end

end