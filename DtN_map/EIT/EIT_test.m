% [p,t,b] = squaregrid(m,n) % create grid of N=mn nodes to be listed in p
% generate mesh of T=2(m-1)(n-1) right triangles in unit square
m=81; n=81; % includes boundary nodes, mesh spacing 1/(m-1) and 1/(n-1)

% computing the grid
[x,y]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1)); % matlab forms x and y lists
p = [x(:),y(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodes

dx = x(2,1) - x(1,1);
%% Building the Triangle List
t = [1,  2,m+2;...
     1,m+2,m+1]; % 3 node numbers for two triangles in first square
t = kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
% now t lists 3 node numbers of 2(m-1) triangles in the first mesh row
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
% final t lists 3 node numbers of all triangles in T by 3 matrix 

%% Volume indices
index_matrix = [1:n*m];
index_matrix = reshape(index_matrix,m,n);

vol_idx = reshape(index_matrix(2:end-1, 2:end-1), (n-2)*(m-2),1) ;

%% Boundary indeces
bdy_idx = [index_matrix(2:end,1);        index_matrix(end,2:end)';...
           index_matrix(end-1:-1:1,end); index_matrix(1,end-1:-1:1)'];

%% proallocating the matrix
% [K,F] = assemble(p,t) % K and F for any mesh of triangles: linear phi's
N=size(p,1);
T=size(t,1); % number of nodes, number of triangles
% p lists x,y coordinates of N nodes, t lists triangles by 3 node numbers

P = 1 + 10*phantom('Modified Shepp-Logan',m);
centroid_T = squeeze(sum(reshape(p(t,:),T,3,2), 2)/3);
% we exchange the positions of x and y given that they were not 
% coming from a meshgrid
sigma = interp2(y,x,P.',centroid_T(:,2),centroid_T(:,1));

[K, ~, ~, ~] = stiff_assembly(m, n, sigma);

%% Compute the data
% [Kb,Fb] = dirichlet(K,F,b) % assembled K was singular! K*ones(N,1)=0
% Implement Dirichlet boundary conditions U(b)=0 at nodes in list b

F = zeros(N,length(bdy_idx));
Dirichlet_data = eye(length(bdy_idx),length(bdy_idx));

F(bdy_idx,:) = Dirichlet_data;

% reduced Stiffness matrix
Kb = K(vol_idx, vol_idx); 

Fb = -K(vol_idx,bdy_idx)*Dirichlet_data; % Stiffness matrix Kb (sparse format) and load vector Fb

% Solving for the vector U will produce U(b)=0 at boundary nodes
U=Kb\Fb;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N

Sol = zeros(N,length(bdy_idx));
Sol(bdy_idx,:) = F(bdy_idx,:);
Sol(vol_idx,:) = U;

Flux = K*Sol/dx;
DtN = Flux(bdy_idx, :); 

%starting point 
Sigma_0 = 1+0*sigma;

%% test withe the function
[l2, grad] = misfit(Sigma_0, DtN, m, n);

J = @(x) misfit(x, DtN, m, n);

options = optimoptions('fminunc','Algorithm','quasi-newton',...
                        'SpecifyObjectiveGradient',true,...
                        'MaxIterations', 10000,...
                        'OptimalityTolerance', 1e-9, ...
                        'Display', 'iter-detailed');
                    
% options = optimoptions('fminunc','Algorithm','quasi-newton',...
%                         'SpecifyObjectiveGradient',true,...
%                         'CheckGradients', true,...
%                         'FiniteDifferenceType', 'central',...
%                         'MaxIterations', 10000,...
%                         'Display', 'iter-detailed');

[x,fval,exitflag,output] = fminunc(J,Sigma_0,options);


%% Assembly of the Mass and differentiation matrices

N = size(p,1);
T = size(t,1);

M = sparse(N, T);
Mass = sparse(N,N); 

% local mass matrix
MK = 1/12*[ 2 1 1;...
            1 2 1;...
            1 1 2];

%% Assembly the matrix
for e=1:T  % integration over one triangular element at a time
  % row of t = node numbers of the 3 corners of triangle e
  nodes = t(e,:);
  
  % 3 by 3 matrix with rows=[1 xcorner ycorner] 
  Pe = [ones(3,1),p(nodes,:)]; 
  % area of triangle e = half of parallelogram area
  Area = abs(det(Pe))/2; 
  % columns of C are coeffs in a+bx+cy to give phi=1,0,0 at nodes
  C=inv(Pe); 
  % now compute 3 by 3 Ke and 3 by 1 Fe for element e
  grad=C(2:3,:);
  % element matrix from slopes b,c in grad
  M_loc=ones(3,1)*Area/3;
  
  M(nodes,e)= M(nodes,e)+M_loc; % add Ke to 9 entries of global K
  Mass(nodes, nodes) = Mass(nodes, nodes) + Area*MK;
  
end   % all T element matrices and vectors now assembled into K and F

sigma_reconstructed = Mass\(M*x);
% 
% % Plotting 
% Plot the FEM approximation U(x,y) with values U_1 to U_N at the nodes 
figure(1); clf();
trisurf(t,p(:,1),p(:,2),0*p(:,1),sigma_reconstructed,'edgecolor','none','facecolor','interp');
view(2),axis equal, colorbar()
% removing ticks and labels
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gcf, 'Position',  [0, 0, 1000, 800])
set(gca,'FontSize', 24)
grid off
axis off
top = 3.22;
bottom = 0.5;
caxis manual

caxis([bottom top]);
% 

% save('reconstructed.dat', 'sigma', 't', 'M', 'p', 'x')


