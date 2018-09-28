%The purpose of this script is to do a fit to of the sum-rate distribution to see how it fits a normal
%distribution by the central limit theorem

rv_len = 10000;
l = 2;  %group size
N = l;  %num TX ants

sigma_n = 0.1;

rho_max = 2;
rho_min = 1;
theta = pi/4;


h_norm_vec = (rho_max-rho_min)*rand(rv_len,l);
h_orth_mat = cos(theta)*(rho_max-rho_min)*rand(rv_len,l,l);

clear c;
for h_norm_idx = 1:l
    intf_sum = zeros(rv_len,1);
    for h_orth_idx = 1:l
        if h_orth_idx == h_norm_idx
            tmp = 0;
        else
            intf_sum = intf_sum + (h_orth_mat(:,h_norm_idx, h_orth_idx).^2);
        end
    end
    sinr(:,h_norm_idx) = h_norm_vec(:,h_norm_idx)./(intf_sum + sigma_n);
    c(:,h_norm_idx) = log2(1+sinr(:,h_norm_idx));
end

if l>1
    sum_rate = sum(c')';
else
    sum_rate = c;
end


[muHat,sigmaHat,muCI,sigmaCI] = normfit(sum_rate)     
