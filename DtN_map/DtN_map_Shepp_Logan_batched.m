L = 1; 
Nx = 2^12; 
dx = L/Nx; 
x0 = [0:dx:L];
[xx,yy] = meshgrid(x0,x0);

nBatches = 128;

%% creating the discretization Grid for the Shepp-Logan phantom
NxSL = 1000;
dxSL = L/NxSL; 
xSL = [0:dxSL:L];
[xxSL,yySL] = meshgrid(xSL,xSL);

P = phantom('Modified Shepp-Logan',NxSL+1);
sigma_fun =  @(x,y) 1 + 0.01*(2+1.8*sin(2*pi*x))./(2+1.8*cos(2*pi*y))+(2+sin(2*pi*y))./(2+1.8*sin(2*pi*x)) + ...
                    10*interp2(xxSL,yySL, P, x,y );

% figure(1);
% sigma_plot = sigma_fun(xxSL, yySL);
% imagesc(xSL,xSL,sigma_plot); axis equal; colorbar();
% set(gca,'FontSize',18)

                
%% ADJUST FREQUENCY HERE
% sigma_fun = @(x,y)(2+1.8*sin(2*pi*x))./(2+1.8*cos(2*pi*y))+(2+sin(2*pi*y))./(2+1.8*sin(2*pi*x));

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

NdirData = length(dirichlet_index);
DirichletData = speye(NdirData);

batchSize = NdirData/nBatches;

flux_d_total = [];

for ii = 1 : nBatches

	Soln = zeros((Nx+1)^2, batchSize);
	Soln(dirichlet_index,:) =  DirichletData(:,...
                               (ii-1)*batchSize+1:ii*batchSize);

    LoadVec = -StiffA(:,dirichlet_index)*DirichletData (:,...
                                (ii-1)*batchSize+1:ii*batchSize);
                           
    LoadVec = LoadVec(unknown_index, :);

    Soln(unknown_index,:)= StiffA(unknown_index,unknown_index)\LoadVec;

    Flux = StiffA*Soln/dx;
    flux_d_total = [flux_d_total Flux(dirichlet_index, :)]; 

end

% save(['test_DE/DtN_',num2str(Nx),'.mat'], 'flux_d_total');
% set(gca,'FontSize',18)
% imagesc(flux_d_total)
% colorbar