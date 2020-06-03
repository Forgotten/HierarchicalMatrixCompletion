L = 1; 
Nx = 2^7; 
dx = L/Nx; 
x0 = [0:dx:L];
[xx,yy] = meshgrid(x0,x0);


%% ADJUST FREQUENCY HERE

sigma_fun = @(x,y)(2+1.8*sin(2*pi*x))./(2+1.8*cos(2*pi*y))+(2+sin(2*pi*y))./(2+1.8*sin(2*pi*x));

StiffA = zeros((Nx+1)^2,(Nx+1)^2);
soln = zeros(Nx+1,Nx+1); soln = soln(:);

% stiff_local_l = [1/2,0,0,-1/2;0,0,0,0;0,0,1/2,-1/2;-1/2,0,-1/2,1];
% stiff_local_r = [1/2,-1/2,0,0;-1/2,1,-1/2,0;0,-1/2,1/2,0;0,0,0,0];
stiff_local_l = [1,-1/2,0,-1/2;-1/2,1/2,0,0;0,0,0,0;-1/2,0,0,1/2];
stiff_local_r = [0,0,0,0;0,1/2,-1/2,0;0,-1/2,1,-1/2;0,0,-1/2,1/2];
sigma_l = sigma_fun(xx(1:end-1,1:end-1) + dx/3,yy(1:end-1,1:end-1) + 2*dx/3);
sigma_r = sigma_fun(xx(1:end-1,1:end-1) + 2*dx/3,yy(1:end-1,1:end-1) + dx/3);

for km = 1:Nx % x direction coordinate
    for kn = 1:Nx % y direction coordinate
        indices = index_map(km,kn,Nx);
        StiffA(indices,indices) = StiffA(indices,indices) + sigma_l(km,kn)*stiff_local_l;
        StiffA(indices,indices) = StiffA(indices,indices) + sigma_r(km,kn)*stiff_local_r;
    end
end

index_matrix = [1:(Nx+1)^2];
index_matrix = reshape(index_matrix,Nx+1,Nx+1);
dirichlet_index = [index_matrix(2:end,1);index_matrix(end,2:end)';...
                   index_matrix(end-1:-1:1,end);index_matrix(1,end-1:-1:1)'];
dirichletdata = zeros(length(dirichlet_index),1);


unknown_index = index_matrix(2:end-1,2:end-1); 
unknown_index = unknown_index(:);


flux_d_total = [];

for km = 1:length(dirichlet_index)

    dirichletdata(km) = 1;
    soln(dirichlet_index) = dirichletdata;
    load_vec = -StiffA(:,dirichlet_index)*soln(dirichlet_index); 
    load_vec = load_vec(unknown_index);

    % figure(3);mesh(reshape(load,Nx,Nx));pause;
    soln(unknown_index)= StiffA(unknown_index,unknown_index)\load_vec;
    flux = StiffA*soln/dx;
    flux_d = flux(dirichlet_index); 
    flux_d_total = [flux_d_total,flux_d];

    soln = reshape(soln,Nx+1,Nx+1);

    %figure(3);mesh(reshape(soln,Nx+1,Nx+1));pause;

    dirichletdata(km) = 0; 
    soln = zeros(Nx+1,Nx+1); soln = soln(:);

end


%save(['test_DE/DtN_',num2str(Nx),'.mat'], 'flux_d_total');



