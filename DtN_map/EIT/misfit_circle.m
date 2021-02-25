function [l2, grad] = misfit_circle(sigma_vec, DtN, t,p, bdy_idx,vol_idx, bdy_triang_idx) 
% This is the misfit for a ciruclar geometry 
if nargin > 6 
    clipping = 1;
else
    clipping = 0;
end
plot = 0;

%% Compute the DtN map for the current conductivity

% assemble the matrix
K = stiff_assembly_from_triang(t, p, sigma_vec);

N = size(p,1);  % number of nodes
T = size(t,1);  % number of triangles

% build the Rhs for the Dtn map
F = zeros(N,length(bdy_idx));
bdy_data = eye(length(bdy_idx),length(bdy_idx));

F(bdy_idx,:) = bdy_data;

% reduced Stiffness matrix
Kb = K(vol_idx, vol_idx); 

Fb = -K(vol_idx,bdy_idx)*bdy_data; % Stiffness matrix Kb (sparse format) and load vector Fb

% TODO solve this by pieces 
% Solving for the vector U will produce U(b)=0 at boundary nodes
U=Kb\Fb;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N

Sol = zeros(N,length(bdy_idx));
Sol(bdy_idx,:) = F(bdy_idx,:);
Sol(vol_idx,:) = U;

Flux = K*Sol;
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
view(2),axis equal,colorbarsype
end

%% Assembly the matrix we needed to extract the Neumann data
Kx = sparse(T,N);
Ky = sparse(T,N);
Surf = sparse(T,T);

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
  Surf(e,e) = Area;
  
end   % all T element matrices and vectors now assembled into K and F


Sol_adj_x = Surf\(Kx*Sol_adj);
Sol_adj_y = Surf\(Ky*Sol_adj); 

Sol_x = Surf\(Kx*Sol);
Sol_y = Surf\(Ky*Sol);

grad = Surf*sum(Sol_adj_x.*Sol_x + Sol_adj_y.*Sol_y, 2);

% Test if this is accurate: We clip the gradient at the boundaty, where 
% we knoe wth exact value of the conductivity
if clipping
    grad(bdy_triang_idx) = 0;
end

if plot
figure(10); clf();
scatter(centroid_T(:,1), centroid_T(:,2),[], grad);
colorbar();
end

l2 = 0.5*sum(sum((DtN_i - DtN).^2));


