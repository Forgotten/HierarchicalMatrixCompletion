function [f] = LTE_cell(epsilon_num,r_number,v_number,freq,incoming)
global Nr Nv epsilon dv dr v0 rr sigma_x

if nargin == 0
	epsilon_num = 0;
    r_number = 50;
    v_number = 60;
    freq = 10;
    incoming = zeros(v_number,1);
end

Lr = 1; dr = 1/r_number; r0 = 0:dr:Lr; Nr = length(r0);
dv = 2/v_number; v0 = -1+dv/2:dv:1-dv/2; Nv = length(v0);
[rr,vv] = meshgrid(r0,v0);

epsilon = 2^(-epsilon_num);

sigma_fun = @(km,x)cos(km*pi*x)+1.5;

sigma_x = sigma_fun(freq,r0);%-sigma_fun(freq,r0)+1;

%% initial 
f = ones(size(rr));
data_pre = zeros(size(f));
for kv = 1:Nv
    v = v0(kv);
    if v>=0
        data_pre(kv,1) = incoming(kv);
%         data_pre(kv,1) = 1;
    else
        data_pre(kv,end) = incoming(kv);
%         data_pre(kv,end) = 0;
    end
end
data_pre = data_pre(:);

f = gmres(@MAB_multi,data_pre,250,1e-9);

% consider the outgoing and incoming
f = reshape(f,Nv,Nr);

f = [f(v0<0,1),f(v0>0,end)];
f = f(:);

% save(['cell/ep_',num2str(epsilon_num),'_freq_',num2str(freq),'.mat']);
% mesh(f);pause;
% 
% save(['test/ep_',num2str(epsilon_num),'_freq_',num2str(freq),'.mat']);

end


function Mb = MAB_multi(p)
global epsilon Nr Nv dr v0 sigma_x
e = ones(Nv,1);
v = v0(:);

p = reshape(p,Nv,Nr);
Rp = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kv = 1:Nv
    p_temp = p(kv,:);
    
    Rp_pre = (p_temp(2:end) - p_temp(1:end-1))/dr;
    
    if v(kv)>=0
        Rp(kv,2:end) = v(kv)*Rp_pre;
        boundary_D(kv,1) = p(kv,1);
    else
        Rp(kv,1:end-1) = v(kv)*Rp_pre;
        boundary_D(kv,end) = p(kv,end);
    end
end

BCp(:,2:Nr-1) = p(:,2:Nr-1) - e*e'*p(:,2:Nr-1)/Nv;

Mb = epsilon*Rp + (e*sigma_x).*BCp + boundary_D;

Mb = Mb(:);

end

