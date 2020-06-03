function Albedo(ep_num)
% denote f_i^- = A_i f_i^+ + B_i f_{i+1}^-
% denote f_{i+1}^+ = C_i f_i^+ + D_i f_{i+1}^-
Nv = 60; Nr = 50;
ep = ep_num; freq = 2;

% left = []; right = [];
f = [];
for km =1:Nv
    km
    incoming = zeros(Nv,1);
    incoming(Nv-km+1)=1;
    f_temp = LTE_cell(ep,Nr,Nv,freq,incoming);
    f = [f,f_temp];
%     left = [left,gg(:,1)];
%     right = [right,gg(:,2)];
%     AB(:,km) = gg(1:Nv/2,1);
%     AB_full(:,km) = gg(:,1);
%     CD(:,km) = gg(Nv/2+1:end,2);
%     CD_full(:,km) = gg(:,2);
end

save(['data/albedo_ep_',num2str(ep),'_freq_',num2str(freq),'.mat']);
% AA = AB(:,1:Nv/2); BB = AB(:,Nv/2+1:end);
% CC = CD(:,1:Nv/2); DD = CD(:,Nv/2+1:end);

% save(['cell/albedo_ep_0_freq_2.mat']);
% save('cell/albedo_ep_0_sigma_1.mat');
end