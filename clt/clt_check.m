%The purpose of this script is to do a fit to of the sum-rate distribution to see how it fits a normal
%distribution by the central limit theorem

rv_len = 10000;
l = 4;  %group size
N = l;  %num TX ants

sigma_n = 0.1;

%The columns and depth dimension represent a given user.
%The rows are indices betwen users.
%Indices are (column,depth,row)
%Note that each user has N independent vectors. There are l different users.
h_mat = (1/sqrt(2*N))*(randn(rv_len,l,l) + i*randn(rv_len,l,l));

clear h_norm_orth_mat;
for h_norm_idx = 1:l
    for h_orth_idx = 1:l
        h_norm_orth_mat(:,h_norm_idx,h_orth_idx) = abs(sum((h_mat(:,:,h_norm_idx).*conj(h_mat(:,:,h_orth_idx)))'));
    end
end

for h_norm_idx = 1:l
    intf_sum = zeros(rv_len,1);
    for h_orth_idx = 1:l
        if h_orth_idx == h_norm_idx
            tmp = 0;
        else
            intf_sum = intf_sum + (h_norm_orth_mat(:,h_norm_idx, h_orth_idx).^2);
        end
    end
    sinr(:,h_norm_idx) = h_norm_orth_mat(:,h_norm_idx, h_norm_idx)./(intf_sum + sigma_n);
    c(:,h_norm_idx) = log2(1+sinr(:,h_norm_idx));
end

sum_rate = sum(c')';


[muHat,sigmaHat,muCI,sigmaCI] = normfit(sum_rate)     
