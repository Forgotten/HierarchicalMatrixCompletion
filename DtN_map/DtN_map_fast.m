% fast version of the DtN map calculation. 
% in this case we generate the sparse matrix in a vectorized form
% in addition we 

L = 1; 
Nx = 2^9; 
dx = L/Nx; 
x0 = [0:dx:L];
[xx,yy] = meshgrid(x0,x0);

%% ADJUST FREQUENCY HERE
sigma_fun = @(x,y)(2+1.8*sin(2*pi*x))./(2+1.8*cos(2*pi*y))+(2+sin(2*pi*y))./(2+1.8*sin(2*pi*x));

%sigma_fun = @(x,y)(10+9.8*sin(9*pi*(x+y)))./(2+1.8*cos(2*pi*y))+(2+sin(2*pi*y))./(2+1.8*sin(2*pi*x))+ ...
    %             5*(exp(-0.5*((x-0.1).^2+ (y-0.3).^2)));


StiffA = sparse((Nx+1)^2,(Nx+1)^2);
soln = zeros(Nx+1,Nx+1); soln = soln(:);

stiff_local_l = [1,-1/2,0,-1/2;-1/2,1/2,0,0;0,0,0,0;-1/2,0,0,1/2];
stiff_local_r = [0,0,0,0;0,1/2,-1/2,0;0,-1/2,1,-1/2;0,0,-1/2,1/2];
sigma_l = sigma_fun(xx(1:end-1,1:end-1) + dx/3,yy(1:end-1,1:end-1) + 2*dx/3);
sigma_r = sigma_fun(xx(1:end-1,1:end-1) + 2*dx/3,yy(1:end-1,1:end-1) + dx/3);

%% Building the Sparse Matrix in vectorized form

[idxJ,idxI] = meshgrid([1:Nx],[1:Nx]);
IdxA = index_mapVec(idxI(:),idxJ(:),Nx);

IdxA2 = repmat(reshape(IdxA.', 1,4*Nx^2), 4, 1);

IdxAJ = reshape(IdxA2, 4, 4, Nx^2);
IdxAI = permute(IdxAJ,[2,1,3]);

% puting them in vector form
IdxJvec = IdxAJ(:);
IdxIvec = IdxAI(:);

% computing the entries and reshaping them 
SigmaL = kron(sigma_l, stiff_local_l);
SigmaR = kron(sigma_r, stiff_local_r);

SigmaL2 = reshape(SigmaL, 4, Nx, 4, Nx);
SigmaLord = permute(SigmaL2, [1,3,2,4]);

SigmaR2 = reshape(SigmaR, 4, Nx, 4, Nx);
SigmaRord = permute(SigmaR2, [1,3,2,4]);

StiffA = sparse(IdxIvec,IdxJvec, SigmaLord(:) + SigmaRord(:), (Nx+1)^2, (Nx+1)^2);

%% getting the different indices for the matrix

index_matrix = [1:(Nx+1)^2];
index_matrix = reshape(index_matrix,Nx+1,Nx+1);

dirichlet_index = [index_matrix(2:end,1);index_matrix(end,2:end)';...
                   index_matrix(end-1:-1:1,end);index_matrix(1,end-1:-1:1)'];

dirichletdata = zeros(length(dirichlet_index),1);

unknown_index = index_matrix(2:end-1,2:end-1); 
unknown_index = unknown_index(:);

%% preparing the RHS in vectorized form
NdirData = length(dirichlet_index);
DirichletData = speye(NdirData);

Soln = zeros((Nx+1)^2, NdirData);
Soln(dirichlet_index,:) =  DirichletData;

LoadVec = -StiffA(:,dirichlet_index)*DirichletData; 
LoadVec = LoadVec(unknown_index, :);

%% solve the system
Soln(unknown_index,:)= StiffA(unknown_index,unknown_index)\LoadVec;

% compute the flux
Flux = StiffA*Soln/dx;
Flux_d_total = Flux(dirichlet_index, :); 

% % save the flux 
% save(['test_DE/DtN_',num2str(Nx),'.mat'], 'flux_d_total');

imagesc(Flux_d_total(10:end/2-10,end/2+10:end-10))
colorbar