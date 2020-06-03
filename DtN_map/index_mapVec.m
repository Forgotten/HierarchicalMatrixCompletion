function indices = index_mapVec(km,kn,Nx_H)
% vectorized version for the creation of the indez
    % we suppose that km and kn are vectors with the same dimension
    assert(size(km,1) == size(kn,1), "indices do not have the same dimension");
    assert(size(km,2) == 1, "indices needs to be a vector");
    
    indices = [(km-1)*(Nx_H+1)+kn, (km-1)*(Nx_H+1)+kn+1, ...
               (km)*(Nx_H+1)+kn+1, (km)*(Nx_H+1)+kn];
           
           
end
