classdef PLR_matrix < handle
    % PLR_matrix : this class contains a constructor and method for
    % Partitioned Low Rank (PLR) matrices using the randomized SVD.
    %
    %
    % A PLR matrix has the following properties:
    % - H_struct : a hierarchical, recursive structure that contains the
    %              size of the current block and its position (top left
    %              entry) in the larger block it is a part of. It also
    %              contains and identifier that says whether the block is
    %              in hierarchical form (has been divided) or in compressed
    %              form. If the block was compressed,  the structure also
    %              contains it's rank, if not the rank is 0. Finally, the
    %              structure contains the data of each block: either the
    %              low rank matrices that approximate it, or another
    %              structure that contains the info about its children.
    % - n        : number of rows of the original matrix
    % - m        : number of columns of the original matrix
    % - r        : maximum rank allowed for the compression
    % - e        : maximum error allowed in compressing an individual
    %              block, in the L2 norm
    %
    % Rosalie Belanger-Rioux <rbr@math.harvard.edu> and Leonardo Zepeda <zepeda@mit.edu>, c. 2014

    
    properties
        H_struct    % struc to store the PLR matrix
        n           % number or rows of the matrix to be compresed
        m           % number of columns of the matric to be compresses
        r           % rank fo the matrix
        e           % error on the approximation 
        norm
        U_sparse    % extra fields to save the sparse factorization
        V_sparse    % extra fields to save the sparse factorization
        form        % in which form the matrix is saved h for H-matrix or s for sparse
        nnz_sparse  % number of non zero elements of the sparse matrix
    end
    
    methods
        %class constructor
        % M_plr = PLR_matrix(M, max rank, error tolerance, norm to compress the matrix)
        function obj = PLR_matrix(D, max_rank, epsilon, normtouse)
            if  nargin == 1   %only one element
                if isa(D,'PLR_matrix') % if that element is a PLR_matrix just copy it
                    obj.H_struct = D.H_struct;
                    obj.n = D.n;
                    obj.m = D.m;
                    obj.r = D.r;
                    obj.e = D.e;
                    obj.norm = D.norm;
                end
                max_rank  = floor(size(D,1)/16);  % a 16th of the full dimension
                epsilon   = 1e-6; %single precision by default
                normtouse ='l2'; % by default the approximation uses the l^2 norm for the compression
            end
            if nargin < 5    % by default we use the induced l^2 norm to decide the compression rate
                normtouse ='l2';
            end
            % we may need to perform a small modification we need to feed
            % the size of the global matrix
            obj.H_struct = dense2H(D, max_rank, epsilon,1,1,normtouse, [size(D,1) size(D,2)]);
            obj.n        = size(D,1);
            obj.m        = size(D,2);
            obj.r        = max_rank;
            obj.e        = epsilon;
            obj.norm     = normtouse;
            [obj.U_sparse, obj.V_sparse] = obj.PLR_2_sparse();
            obj.nnz_sparse = nnz(obj.U_sparse) + nnz(obj.V_sparse);
            obj.form     = 's';   
        end % class constructor
        
        % container for the recursive function
        function [hist_size szs nums mydisp]= histogram(obj, flag, opt, name)
            % function [comp, hist_size, szs, nums, mydisp]= histogram(obj, flag, opt, plt, name)
            % Compute important information about the compressed matrix
            %
            % Inputs:
            % - obj  : a PLR_matrix_rand
            % - flag : string, which kind of plot to create:
            %      -- 'hist'  : use the built-in Matlab histogram (broken?)
            %      -- 'graph' : sort the square roots of nb of points for
            %                   each block, plot those
            %      -- 'draw'  : draw the outlines of each block inside the
            %                   matrix, and color each block by its
            %                   compressed rank
            %      -- 'freq'  : plot how many blocks of each size there are
            % - opt  : if we 'draw', then opt='log' means a log color
            %          scale, any other opt or no opt is the standard
            %          scale. If flag='freq' instead, then opt toggles thep
            %          lotting of the 'ideal' frequencies: opt=1 or no opt
            %          means 'ideal' frequencies for diagonal blocks. An
            %          opt of 2 means 'ideal' for blocks jsut off the
            %          diagonal, any other opt does not plot ideal
            %          frequencies.
            % - plt  : if flag='freq', plt is the color for plotting the
            %          actual frequencies (ideals are plotted blue for the
            %          diagonal, black for off the diagonal), if no plt
            %          then default is magenta ('m'). Red is 'r', green
            %          'g', yellow 'y', cyan 'c'.
            % - name : string put in the legend to describe the actual
            %          frequencies, default is 'block size frequencies'.
            %
            % Outputs:
            % - comp           : approximate complexity of multiplying this
            %                    compressed matrix with a vector: 4 times
            %                    the sum over all blocks of the product of
            %                    the rank with the number of rows
            % - hist_size(1,:) : number of rows of each block
            % - hist_size(2,:) : number of colns of each block
            % - hist_size(3,:) : rank of each block
            % - hist_size(4,:) : top-most row index of each block
            % - hist_size(5,:) : left-most column index of each block
            % - szs            : contains the unique products of sizes of
            %                    blocks
            % - nums           : number of blocks of each size in szs
            % - mydisp         : matrix containing a number corresponding
            %                    to either the rank of the leaf block that
            %                    point is in, or 1/2 if the point is on the
            %                    edge of a block. That matrix is used when
            %                    we 'draw' the outlines of the blocks with
            %                    their ranks
            %

            mydisp    = zeros(obj.n,obj.m);
            N         = min(obj.n,obj.m);
            hist_size = histogram_size(obj.H_struct);
            szs       = unique(hist_size(1,:).*hist_size(2,:));
            nums      = zeros(size(szs));
            for ins = 1:length(szs)
                nums(ins) = nnz(hist_size(1,:).*hist_size(2,:)==szs(ins));
            end
            if nargin > 1
                if strcmp(flag, 'hist') % we show the built-in histograme of the sizes
                    hist(hist_size);
                elseif strcmp(flag, 'graph');  % we show just a graph of the squared root of the block sizes
                    plot(sort(sqrt(hist_size)));
                elseif strcmp(flag, 'draw'); % draw the blocks with their ranks
                    clf
                    for myb=1:size(hist_size,2)
                        % put in the rank
                        inr = hist_size(4,myb):hist_size(4,myb)+hist_size(1,myb)-1;
                        inc = hist_size(5,myb):hist_size(5,myb)+hist_size(2,myb)-1;
                        mydisp(inr,inc) = hist_size(3,myb);
                        % trace outline
                        mydisp(inr,inc(1))   = 1/2;
                        mydisp(inr,inc(end)) = 1/2;
                        mydisp(inr(1),inc)   = 1/2;
                        mydisp(inr(end),inc) = 1/2;
                    end
                    if nargin>2 && strcmp(opt,'log'),
                        imagesc(log10(mydisp));
                    else imagesc(mydisp);
                    end
                    colorbar;
                    if obj.n == obj.m, axis square;end
                elseif strcmp(flag, 'freq'); % get freq of each size
                    if nargin > 2
                        pltstr=['-' opt 'o'];
                    else
                        pltstr='-bo';
                    end
                    if nargin > 3
                        lgd = name; 
                    else 
                        lgd ='block size frequencies';
                    end
                    
                    L = log2(N/obj.r);
                    loglog((N./2.^(1:L)).^2,[2.^(1:L-1) 2^(L+1)],'-k','LineWidth',2);hold on;
                    %loglog(2.^(2*[1:log2(N/2)]),N./sqrt(2.^(2*[1:log2(N/2)])),'-k');hold on;
                    
                    loglog(szs,nums,pltstr);
                    legend('ideal frequencies per block size',lgd);     grid on;
                    xlabel('block size in number of entries');           ylabel('frequency');
                end
            end
        end  % plot the histogram of the size of the blocks
        
        function [U, V] = PLR_2_sparse(self)
            [U,V] = H2sparse(self.H_struct);
        end
        
        
        function r = mtimes(H,x)
            % mtimes overload function
            if isvector('x')
                if H.form == 'h'
                r = H_mat_vec(H.H_struct,x);
                %r = r(:); %RBR: bad if multiple column vectors
                elseif H.form == 's'
                r = H.U_sparse*(H.V_sparse*x);   % using the sparse factorization form 
                end
            else
                fprintf('M*x, x needs to be a vector');
            end
        end % mat_vec multiplication
        
        function M = PLR_2_dense(obj)
            % wrapper of the recursive function to decompress the matrix
            M = H2dense(obj.H_struct);
        end % PLR to dense matrix decompression
        
    end
    
end

% auxiliray function wrapped inside the class
function M = H2dense(H_struct)
    % function to decompress a PLR or H matrix
    if H_struct.id == 'c' % it's a compressed matrix ( a leaf)
        M = H_struct.Data{1}*H_struct.Data{2};
    else                  % it's a branch so we proceed by recursivity
        M_11 = H2dense(H_struct.Data{1,1});
        M_12 = H2dense(H_struct.Data{1,2});
        M_21 = H2dense(H_struct.Data{2,1});
        M_22 = H2dense(H_struct.Data{2,2});
        M = [M_11 M_12 ;
             M_21 M_22];
    end
    
end



function H = dense2H(M, max_rank, epsilon , starti , startj , normtouse, global_size)
%function to compute an adaptive H-matrix from a dense matrix.
% Inputs : 
% - M        : dense matrix or block
% - max_rank : maximum rank of compressed blocks
% - epsilon  : minimum accuracy per block (L2 norm)
%
% output:
% - H        : structure containing the hierarchical info
%
% in some cases the SVD does not converge, to circumvent this issue we
% handle the exception and if the svd did not converge we split the matrix

    
%     n    = size(M,1); %size of the matrix
%     m    = size(M,2); %size of the matrix
%     
%     error_flag = 0;
% 
%     % try-catch error to handle the possible non convergence of the QR
%     % iteration
%     try
%         [ U, S, V] = svd(M, 0);
%     catch ME
%         ME.stack;
%         error_flag = 1;
%         %fprintf('svd did not converge, trying again \n');
%     end
%     %[ U, S, V, error_flag] = svds(M, min(n,max_rank+1)); % in real life we need to compute just a couple of mat_vec operations
%     %
%     if error_flag == 1 % if didn't converge, try again
%         try 
%             [ U, S, V, error_flag] = svds(M, min(n,max_rank+1)); % in real life we need to compute just a couple of mat_vec operations
%         catch ME 
%             ME.stack;
%             error_flag = 1;
%             fprintf('svd did not converge, we are just dividing the matrix \n');
%         end
%     end
%     

% start randomized SVD 

n     = size(M,1); %size of the matrix
m     = size(M,2); %size of the matrix
nRand = min( max_rank + 5, min(n,m) ); % number of random vectors to use
MatRand = randn(nRand, n);
R = MatRand*M;

error_flag = 0;

try
    [Q, ~, ~] = svd(R', 0);
catch ME
    ME.message
    ME.stack
    error_flag = 1;
end

Q = Q(:, 1:min(max_rank+1,min(n,m))); % we use one more than max_rank, so we can tell with the singular values whether we've converged
T = M*Q;

try
    [U, S, W] = svd(T, 0);
catch ME
    ME.message
    ME.stack
    error_flag = error_flag + 2;
end

V = Q*W;

% finish randomized SVD
   
    ind_size = min(size(S)); % real size of S have to use min because the matrix may not be square
    s = diag(S);
    
    % we define the cut-off on the singular values depending on the norm
    % used
    if strcmp(normtouse, 'l2') % L2 norm
        ind_rank = find((epsilon-s)>0, 1 ) - 1; % empty minus 1 = empty
        if isempty(ind_rank)
            ind_rank =ind_size;
        end
    elseif strcmp(normtouse, 'fro')  % Frobenius norm
        lowsing = find(cumsum(flipud(s(:)).^2)>epsilon^2,1); % compute the cumulative sum of the eigenvalues
        if isempty(lowsing)  % if the matrix can not be compressed
            lowsing = ind_size;
        end
        ind_rank = ind_size+1-lowsing; %RBR: for frobenius norm, if empty then can throw out al lthe matrix! or ind_rank=1
    else % rise error
        fprintf('The metric for the approximation must be either l2 or fro \n By default it is the l2 norm \n');
    end
    % once the cut-off is selected we save the compressed form, of we make
    % a recursive call / The way the matrix is saved is suboptimal when the
    % matrix is small with high rank, however, it's easy to implement
    if (n <= max_rank || m <= max_rank) &&  error_flag == 0 %RBR: if matrix/block is rectangular
        %Data = {U(:,1:ind_rank)*S(1:ind_rank ,1:ind_rank ), V(:,1:ind_rank )'};
        if n < m
            Data = {speye(n),M};
        else 
            Data = {M, speye(m)};
        end
        rnk  = min(n,m);
        %rnk  = ind_rank;
        id   = 'c';  % we have a small enough matrix, we just need to store it
    elseif   error_flag == 0 && ind_rank<=max_rank
        Data = {U(:,1:ind_rank )*S(1:ind_rank ,1:ind_rank ), V(:,1:ind_rank )'};
        rnk  = ind_rank;
        id   = 'c'; % compressed matrix just need to be stored
    else  % recursive call 
        Data = cell(2,2);
        for ii = 1:2
            if ii == 1
                index_rows = 1:floor(size(M,1)/2);
            else
                index_rows = floor(size(M,1)/2)+1:size(M,1);
            end
            for jj = 1:2
                if jj == 1
                    index_columns = 1:floor(size(M,2)/2);
                else
                    index_columns = floor(size(M,2)/2)+1:size(M,2);
                end
                % splitting the matrix and computing sub matrices low rank decompositions
                % we add the information about the relative position of the sub-matrices. 
                % ( this is needed for the plotting routine)
                sti = starti - 1 + index_rows(1);
                stj = startj - 1 + index_columns(1);  % top left corner of block
                % saving the results of the recursive call within the data structure
                Data{ii,jj} = dense2H(M(index_rows,index_columns), max_rank, epsilon , sti , stj , normtouse,  global_size );
            end
        end
        id  = 'h'; % H-matrix, has to be split and compute lower rank approximations
        rnk = 0;
    end
    H = struct('n',n,'m',m,'id', id, 'Data', {Data},'rnk',rnk,'sti',starti,'stj',startj, 'global_size', global_size);
end % converter from dense to H_matrix struct

function [U,V] = H2sparse(H) % this is the recursive function

    if H.id == 'h' % this is a H-matrix has to call recursively the children
        % recursive call the the other matrices
        [U_1, V_1 ] = H2sparse(H.Data{1,1}); 
        [U_2, V_2 ] = H2sparse(H.Data{1,2});
        [U_3, V_3 ] = H2sparse(H.Data{2,1});
        [U_4, V_4 ] = H2sparse(H.Data{2,2});
        
        U = [U_1 U_2 U_3 U_4];
        V = [V_1;
             V_2;
             V_3;
             V_4;];
    elseif H.id == 'c'
        N = H.global_size(1);
        M = H.global_size(2);
        U = sparse(N,H.rnk);
        V = sparse(H.rnk, M);
        
        U(H.sti:H.sti+H.n-1,:) = H.Data{1} ;
        V(:,H.stj:H.stj+H.m-1) = H.Data{2} ;
    end
    
end


function histo = histogram_size(H)
    %recursive function to obtain the size of all the blocks
    if H.id == 'c'
        histo = [H.n;H.m;H.rnk;H.sti;H.stj];
    elseif H.id == 'h'
        histo1 =  histogram_size(H.Data{1,1});
        histo2 =  histogram_size(H.Data{1,2});
        histo3 =  histogram_size(H.Data{2,1});
        histo4 =  histogram_size(H.Data{2,2});
        histo  = [histo1 histo2 histo3 histo4];
    end
end


function y = H_mat_vec( H , x )
    % Function to compute the mat_vec product between a dense vector and an
    % H-matrix.
    % input : H = H-matrix to be applied
    %         x = dense vector to be multiplied
    % output:
    %         y = dense vector H*x
    
    if(numel(x) == H.m)
        if(H.id == 'c') % matrix is compressed just need to evaluate
            % y = (U*S)*(V.'*x)
            y = H.Data{1}*(H.Data{2}*x);
            
        elseif H.id == 'h' % H matrix, then it's split in smaller matrices.
            y = zeros(H.n,size(x,2)); % RBR: for multiple column vectors
            %        y = zeros(H.n,1); % allocating the answer in memory
            for ii = 1:2
                if ii == 1  % vertical blocks
                    index_rows = 1 : floor(H.n/2);
                else
                    index_rows = floor(H.n/2)+1 : H.n;
                end
                for jj = 1:2   %horizontal blocks
                    if jj == 1
                        index_columns = 1 : floor(H.m/2);
                    else
                        index_columns = floor(H.m/2)+1 : H.m;
                    end
                    % y = [H_11*x_1 + H_12*x_2;
                    %      H_21*x_1 + H_22*x_2 ];
                    y(index_rows,:) = y(index_rows,:) + H_mat_vec(H.Data{ii,jj}, x(index_columns,:) ); % RBR: for multiple column vectors
                    %                y(index_rows) = y(index_rows) + H_mat_vec(H.Data{ii,jj}, x(index_columns) );
                end
            end
        end
        
    else
        fprintf('Dimension mismatch');
        y = 0;
    end
end % H matrix - vector multiplication