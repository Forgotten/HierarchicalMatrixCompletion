function [block, mask, err] = block_complete_o(A, p, min_dim, r, tol)
% function to recursively compress an upper diagonal block
if min_dim >= size(A,1)
    
    block = A;
    mask = ones(size(A));
    err = 0.0;
else
    % we run the optimization several times to check if we reach 
    % the given tolerance
    for ii = 1:10
        [block, mask, err] = block_completion(A, p);
        if err < tol
            break
        end
    end
    % this is the baseline
    for jj = 1:10
        p_temp = p*(0.93)^(kk-1);
        for ii = 1:10
            [block_temp, mask_temp, err_temp] = block_completion(A, p_temp);
            if err_temp < tol
                break
            end
        end
        
        % check if the reconstruction is still above the given tolerance
        if err_temp < tol && nnz(mask_temp) < nnz(mask)
            block = block_temp;
            mask = mask_temp;
            err = err_temp;
        end
    end
end


end