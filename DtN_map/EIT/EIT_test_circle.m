% [p,t,b] = squaregrid(m,n) % create grid of N=mn nodes to be listed in p
% generate mesh of T=2(m-1)(n-1) right triangles in unit square

% add distmesh to the path
addpath('../distmesh');

h = 0.0246;
 h = 0.01;
% computing the mesh
fd = @(p) sqrt(sum(bsxfun(@minus,p,[0, 0.0] ).^2,2)) - 1;

[p,t] = distmesh2d(fd,@huniform,h,[-1,-1;1,1],[]);

bdy_idx = find(abs(fd(p)) < 1e-5);
vol_idx = find(abs(fd(p)) > 1e-5);

% we need to order the dby indixes
bdy_x_pos = find(p(bdy_idx, 1) > 0);
bdy_x_neg = find(p(bdy_idx, 1) <= 0);

[~, sort_pos_idx] = sort(p(bdy_idx(bdy_x_pos),2));
[~, sort_neg_idx] = sort(p(bdy_idx(bdy_x_neg),2));

bdy_idx = [ bdy_idx(bdy_x_pos(flipud(sort_pos_idx)));...
            bdy_idx(bdy_x_neg(sort_neg_idx))];

% we add the indices of the triangles to clip de gradient there
bdy_triang_idx = find(sum(ismember(t, bdy_idx), 2));

N = size(p,1);
T = size(t,1); % number of nodes, number of triangles

centroid = squeeze(sum(reshape(p(t,:),T,3,2), 2)/3);

sigma = @(p) 1 + 2*(sqrt(sum(bsxfun(@minus,p,[   0.1,  0.4] ).^2,2))< 0.3) ...
               + 4*(sqrt(sum(bsxfun(@minus,p,[-0.3, -0.3] ).^2,2))< 0.2);

sigma_vec = sigma(centroid);

figure(1); clf();
trisurf(t,p(:,1),p(:,2),0*p(:,1),sigma(p),'edgecolor','none','facecolor','interp');
view(2),axis equal,colorbar
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
axis off
set(gca,'FontSize',24)

%% Compute the DtN map for the current conductivity

% assemble the matrix
K= stiff_assembly_from_triang(t, p, sigma_vec);

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

Flux = K*Sol;
DtN_i = Flux(bdy_idx, :); 
%starting point 
Sigma_0 = 1+0*sigma_vec;


%% test withe the function
[l2, grad] = misfit_circle(sigma_vec, DtN_i, t,p, bdy_idx,vol_idx, bdy_triang_idx);

J = @(x) misfit_circle(x, DtN_i, t,p, bdy_idx,vol_idx, bdy_triang_idx) ;

options = optimoptions('fminunc','Algorithm','quasi-newton',...
                        'SpecifyObjectiveGradient',true,...
                        'MaxIterations', 10000,...
                        'OptimalityTolerance', 1e-9, ...
                        'Display', 'iter-detailed');
        
% % without clipping to check the gradients

% J = @(x) misfit_circle(x, DtN_i, t,p, bdy_idx,vol_idx);
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

% % Plotting 
% Plot the FEM approximation U(x,y) with values U_1 to U_N at the nodes 
% figure(2); clf();
% trisurf(t,p(:,1),p(:,2),0*p(:,1),M*sigma/Area,'edgecolor','k','facecolor','interp');
% view(2),axis equal,colorbar
% 

% save('sigma_reconstructed_normal.mat', 'sigma_reconstructed', 'sigma', ...
%      't', 'M', 'p', 'x')
% 

