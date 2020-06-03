function LTE_2D_ep_direction(ep_index,para_index)
% This function solves RTE with para_sigma_x and boundary 
% condition data_pre
global Nr Nt epsilon dr theta0 para_sigma_x

L = 0.6; dr = 0.05; x = [0:dr:L];
dtheta = 2*pi/12; theta0 = [-pi+dtheta/2:dtheta:pi-dtheta/2];
epsilon = 1/2^(ep_index); Nr = length(x); Nt = length(theta0);

[xx0,yy0] = meshgrid(x,x); [xx,yy,vv] = meshgrid(x,x,theta0);
para_sigma_x = ones(size(xx0));
if para_index == 2
    circle = (xx0-0.3).^2+(yy0-0.3).^2;
    para_sigma_x = ones(size(xx0)) + 0.1*double(circle<0.04);
end
if para_index == 3
    circle = (xx0-0.3).^2+(yy0-0.3).^2;
    para_sigma_x = ones(size(xx0)) + 0.05*double(circle<0.04);
end
if para_index == 4
    circle = (xx0-0.3).^2+(yy0-0.3).^2;
    para_sigma_x = ones(size(xx0)) + 0.025*double(circle<0.04);
end
if para_index == 5
    circle = (xx0-0.3).^2+(yy0-0.3).^2;
    para_sigma_x = ones(size(xx0)) + 0.0125*double(circle<0.04);
end


boundary_index = zeros(size(xx));
boundary_index(2:end-1,(x==0),(cos(theta0)>0)) = ones(Nr-2,1,Nt/2);
boundary_index(2:end-1,(x==L),(cos(theta0)<0)) = ones(Nr-2,1,Nt/2);
boundary_index((x==L),2:end-1,(sin(theta0)<0)) = ones(1,Nr-2,Nt/2);
boundary_index((x==0),2:end-1,(sin(theta0)>0)) = ones(1,Nr-2,Nt/2);
boundary_index = boundary_index(:);

out_index = zeros(size(xx));
out_index(2:end-1,(x==0),(cos(theta0)<0)) = ones(Nr-2,1,Nt/2);
out_index(2:end-1,(x==L),(cos(theta0)>0)) = ones(Nr-2,1,Nt/2);
out_index((x==L),2:end-1,(sin(theta0)>0)) = ones(1,Nr-2,Nt/2);
out_index((x==0),2:end-1,(sin(theta0)<0)) = ones(1,Nr-2,Nt/2);
out_index = out_index(:); k_out = find(out_index);
% 
% boundary_index = zeros(size(xx));
% boundary_index(1,2:end-1,(cos(theta0)>0)) = ones(1,Nr-2,Nt/2);
% boundary_index(end,2:end-1,(cos(theta0)<0)) = ones(1,Nr-2,Nt/2);
% boundary_index(2:end-1,1,(sin(theta0)>0)) = ones(Nr-2,1,Nt/2);
% boundary_index(2:end-1,end,(sin(theta0)<0)) = ones(Nr-2,1,Nt/2);
% boundary_index = boundary_index(:);
% 
% out_index = zeros(size(xx));
% out_index(1,2:end-1,(cos(theta0)<0)) = ones(1,Nr-2,Nt/2);
% out_index(end,2:end-1,(cos(theta0)>0)) = ones(1,Nr-2,Nt/2);
% out_index(2:end-1,1,(sin(theta0)<0)) = ones(Nr-2,1,Nt/2);
% out_index(2:end-1,end,(sin(theta0)>0)) = ones(Nr-2,1,Nt/2);
% out_index = out_index(:); k_out = find(out_index);

k = find(boundary_index);

f_total = [];
data_pre = zeros(size(boundary_index));
for km = 1:length(k);
    data_pre(k(km)) = 1;
    f = gmres(@MAB_multi,data_pre,150,1e-12);
    f_total = [f_total,f(k_out)];
    f = reshape(f,Nr,Nr,Nt);
    data_pre(k(km)) = 0;
end

save(['data_direction_dense/green_Para',num2str(para_index),'_ep_',num2str(ep_index),'.mat']);
end

function Mb = MAB_multi(p)
global epsilon Nr Nt dr theta0 para_sigma_x
e = ones(Nt,1);

p = reshape(p,Nr,Nr,Nt);
Rp_x = zeros(size(p));
Rp_y = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kt = 1:Nt
    p_temp = p(:,:,kt);
    
   Rp_pre_x = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
   Rp_pre_y = (p_temp(2:end,:) - p_temp(1:end-1,:))/dr;
   
    sinT = sin(theta0(kt));
    cosT = cos(theta0(kt));
    if sinT<=0
        Rp_y(1:end-1,:,kt) = sinT*Rp_pre_y;
         Rp_y(:,1,kt) = zeros(Nr,1);
        Rp_y(:,end,kt) = zeros(Nr,1);
        boundary_D(end,:,kt) = p_temp(end,:);
    else
        Rp_y(2:end,:,kt) = sinT*Rp_pre_y;
        Rp_y(:,1,kt) = zeros(Nr,1);
        Rp_y(:,end,kt) = zeros(Nr,1);
        boundary_D(1,:,kt) = p_temp(1,:);
    end
    if cosT>=0
        Rp_x(:,2:end,kt) = cosT*Rp_pre_x;
       Rp_x(1,:,kt) =  zeros(1,Nr);
        Rp_x(end,:,kt) = zeros(1,Nr);
        boundary_D(:,1,kt) = p_temp(:,1);
    else
        Rp_x(:,1:end-1,kt) = cosT*Rp_pre_x;
        Rp_x(1,:,kt) =  zeros(1,Nr);
        Rp_x(end,:,kt) = zeros(1,Nr);
        boundary_D(:,end,kt) = p_temp(:,end);
    end
end

for kx = 2:Nr
    for ky = 2:Nr
        p_temp = p(kx,ky,:); p_temp = p_temp(:);
        p_temp = p_temp - e*e'*p_temp/Nt;
        BCp(kx,ky,:) = para_sigma_x(kx,ky)*p_temp;
    end
end

Mb = epsilon*Rp_x + epsilon* Rp_y + BCp + boundary_D;

Mb = Mb(:);

end

