function [l2, grad] = misfit(sigma_vec, DtN, m, n) 
clipping = 0;
plot = 0;
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


%% Compute the DtN map for the current conductivity

% assemble the matrix
[K, ~, ~, ~] = stiff_assembly(m, n, sigma_vec);

% build the Rhs for the Dtn map
F = zeros(N,length(bdy_idx));
bdy_data = eye(length(bdy_idx),length(bdy_idx));

F(bdy_idx,:) = bdy_data;

% reduced Stiffness matrix
Kb = K(vol_idx, vol_idx); 

Fb = -K(vol_idx,bdy_idx)*bdy_data; % Stiffness matrix Kb (sparse format) and load vector Fb

% Solving for the vector U will produce U(b)=0 at boundary nodes
U=Kb\Fb;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N

Sol = zeros(N,length(bdy_idx));
Sol(bdy_idx,:) = F(bdy_idx,:);
Sol(vol_idx,:) = U;

Flux = K*Sol/dx;
DtN_i = Flux(bdy_idx, :); 

%% Plotting 
if plot
% Plot the FEM approximation U(x,y) with values U_1 to U_N at the nodes 
figure(3); clf();
trisurf(t,p(:,1),p(:,2),0*p(:,1),Sol(:,2),'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar

figure(4); clf();
trisurf(t,p(:,1),p(:,2),0*p(:,1),Sigma_0,'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar

end
%% Compute the adjoint

Adj_rhs = -(DtN - DtN_i);
F_adj = zeros(N,length(bdy_idx));

F_adj(bdy_idx,:) = Adj_rhs;

Fb_adj = -K(vol_idx,bdy_idx)*Adj_rhs; % Stiffness matrix Kb (sparse format) and load vector Fb

% Solving for the vector U will produce U(b)=0 at boundary nodes
U_adj=Kb\Fb_adj;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N

Sol_adj = zeros(N,length(bdy_idx));
Sol_adj(bdy_idx,:) = F_adj(bdy_idx,:);
Sol_adj(vol_idx,:) = U_adj;

%% Plotting 
% Plot the FEM approximation U(x,y) with values U_1 to U_N at the nodes 
if plot
figure(5); clf();
trisurf(t,p(:,1),p(:,2),0*p(:,1),Sol(:,2),'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar
end

%% Assembly the matrixwe needed to extract the Neumann data
Kx = sparse(T,N);
Ky = sparse(T,N);

for e=1:T  % integration over one triangular element at a time
  % row of t = node numbers of the 3 corners of triangle e
  nodes = t(e,:);
  
  % 3 by 3 matrix with rows=[1 xcorner ycorner] 
  Pe = [ones(3,1),p(nodes,:)]; 
  % area of triangle e = half of parallelogram area
  Area = abs(det(Pe))/2; 
  % columns of C are coeffs in a+bx+cy to give phi=1,0,0 at nodes
  C = inv(Pe); 
  % now compute 3 by 3 Ke and 3 by 1 Fe for element e
  grad=C(2:3,:);
  % element matrix from slopes b,c in grad
  Kx_loc = grad(1,:)*Area;
  Ky_loc = grad(2,:)*Area;
  
  Kx(e,nodes) = Kx(e,nodes)+Kx_loc; % add Ke to 9 entries of global K
  Ky(e,nodes) = Ky(e,nodes)+Ky_loc;
  
end   % all T element matrices and vectors now assembled into K and F


Sol_adj_x = Kx*Sol_adj/Area;
Sol_adj_y = Ky*Sol_adj/Area; 

Sol_x = Kx*Sol/Area;
Sol_y = Ky*Sol/Area;

grad = sum(Sol_adj_x.*Sol_x + Sol_adj_y.*Sol_y, 2)*Area;

% Test if this is accurate: We clip the gradient at the boundaty, where 
% we knoe wth exact value of the conductivity
if clipping
    idxT = [1:T];
    c = reshape(idxT, m-1,2*(n-1));
    bdy_idxT = [idxT(2:end,1);        idxT(end,2:end)';...
                idxT(end-1:-1:1,end); idxT(1,end-1:-1:1)'];
    grad(bdy_idxT) = 0;
end

if plot
figure(10); clf();
scatter(centroid_T(:,1), centroid_T(:,2),[], grad);
colorbar();
end

l2 = 0.5*sum((DtN_i - DtN).^2*dx, [1,2]);


