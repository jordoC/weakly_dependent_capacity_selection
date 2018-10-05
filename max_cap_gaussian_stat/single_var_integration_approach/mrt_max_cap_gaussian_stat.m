%Statistical analysis for maximum sum rate for a collection of users in a fading channel with MRT beamforming.
clear;
rv_len = 1000;
%NOTE: l cannot be = 1!!!
% the sums used below will reduce the random vectors to scalars if l = 1
% (ie. the rows will be summed instead of the columns)
l = 4;  %group size
try
    assert(l>1);
catch me
    error(['Group size l (%i) out of range (l must be > 1). Error details: %s'],...
        l, me.message)
end
N = 4;          %num TX ants
tx_power = 1;   %total tx power in Watts
n = 10;          %size of candidate set

sigma_n = 0.1;  %noise variance


%cell array that contains a circularly symmetric normal rv that represents
% the channel for a given user and given antenna. each user has a vector
% of N RVs to describe the MIMO channel (one for each antenna)
h_arr = cell(n,1);
for arr_idx = 1:n
    h_arr{arr_idx,1} = (1/sqrt(2*N))*(randn(rv_len,N,1) + i*randn(rv_len,N,1));
end

%rows of this variable ar the sus groups, each cell contains random vector for a given user
% that corresponds to its chanel vector.
sus_groups = nchoosek(h_arr,l);

num_groups = nchoosek(n,l);
sus_groups_norm_orth = cell(num_groups,1);


%form the norms and inner products for each user
clear h_norm_orth_mat;
for sus_group_idx = 1:num_groups
    %retuns a column cell array 
    %each column contains a random vector for each user (implemented as an array here)
    sus_group = sus_groups(sus_group_idx,:)';
    for h_norm_idx = 1:l
        for h_orth_idx = 1:l
            %this matrix contains norms and inner products of random channel vectors.
            % the norms are accross the diagonal, while the combinations of inner products
            % between channel vectors are off-diagonal. Matrix is symmetric.
            h_norm_orth_mat(:,h_norm_idx,h_orth_idx) = ...
                abs(sum(((sus_group{h_norm_idx}).*conj(sus_group{h_orth_idx}))')');
        end
    end
    sus_groups_norm_orth{sus_group_idx} = h_norm_orth_mat;
end

%form the random variables for sum rate for each sus group
sus_groups_cap = cell(num_groups,1);
sus_groups_sinr = cell(num_groups,1);
sus_groups_sumrt = cell(num_groups,1);
clear h_norm_orth_mat;
clear cap;
clear sinr;
for sus_group_idx = 1:num_groups
    %this matrix is for RVs representing norm and inner products between channel vectors
    h_norm_orth_mat = sus_groups_norm_orth{sus_group_idx};
    for h_norm_idx = 1:l
        intf_sum = zeros(rv_len,1);
        for h_orth_idx = 1:l
            if h_orth_idx == h_norm_idx
                tmp = 0;
            else
                intf_sum = intf_sum + (h_norm_orth_mat(:,h_norm_idx, h_orth_idx).^2);
            end
        end
        sinr(:,h_norm_idx) = (tx_power/N)*h_norm_orth_mat(:,h_norm_idx, h_norm_idx)./...
            ((tx_power/N)*intf_sum + sigma_n);
        cap(:,h_norm_idx) = log2(1+sinr(:,h_norm_idx));
    end
    sus_groups_sinr{sus_group_idx} = sinr;
    sus_groups_cap{sus_group_idx} = cap;
    sus_groups_sumrt{sus_group_idx} = sum(cap')';
end

cov_mat = zeros(num_groups);
corr_mat = zeros(num_groups);
mean_vec = zeros(num_groups,1);
for sus_group_idx_1 = 1:num_groups
    mean_vec(sus_group_idx_1,1) = mean(sus_groups_sumrt{sus_group_idx_1,1});
    for sus_group_idx_2 = sus_group_idx_1+1:num_groups
        pairwise_cov_mat = cov(sus_groups_sumrt{sus_group_idx_1,1}, ...
            sus_groups_sumrt{sus_group_idx_2,1});
        cov_mat(sus_group_idx_1,sus_group_idx_1) = pairwise_cov_mat(1,1);
        cov_mat(sus_group_idx_2,sus_group_idx_2) = pairwise_cov_mat(2,2);
        cov_mat(sus_group_idx_1,sus_group_idx_2) = pairwise_cov_mat(1,2);
        cov_mat(sus_group_idx_2,sus_group_idx_1) = pairwise_cov_mat(2,1);

        corr_mat_denom = sqrt(cov_mat(sus_group_idx_1,sus_group_idx_1))*...
            sqrt(cov_mat(sus_group_idx_2,sus_group_idx_2));
        corr_mat(sus_group_idx_1,sus_group_idx_1) = 1;
        corr_mat(sus_group_idx_2,sus_group_idx_2) = 1;
        corr_mat(sus_group_idx_1,sus_group_idx_2) = cov_mat(sus_group_idx_1,sus_group_idx_2)/...
            corr_mat_denom;
        corr_mat(sus_group_idx_2,sus_group_idx_1) = cov_mat(sus_group_idx_2,sus_group_idx_1)/...
            corr_mat_denom;
    end
end



%model the sum rate as a normal distribution by central limit theorem.
% we use the properties of dependent random variables to come up with the
% pdf for the maximum sum rate accross all the sus groups.

%build up to evaluate the multi-var cdf, marginal pdf on
mu = mean_vec(1);
sigma = sqrt(mean(diag(cov_mat)));

rho_vec = [];
for off_diag_idx = 1:num_groups
    rho_vec = [rho_vec' diag(corr_mat,off_diag_idx)']';
end

rho = mean(rho_vec);

f_n = @(z) num_groups*(normcdf(z,0,1).^(num_groups-1)).*normpdf(z,0,1);
F_n = @(z) normcdf(z,0,1).^num_groups;

g_int_arg = @(x,z) (1/(sigma*sqrt(1-rho)))*f_n((((x-mu)/sigma) + sqrt(rho)*z)/(sqrt(1-rho))).*normpdf(z,0,1);
G_int_arg = @(x,z) (F_n((((x-mu)/sigma) + sqrt(rho)*z)/(sqrt(1-rho)))).*normpdf(z,0,1);


clear g_n;
clear G_n;
x = mu+(-2*sigma:sigma/30:6*sigma); %bounds of grid are estimated from first rv
g_n(:,1) = x;
G_n(:,1) = x;
for i = 1:length(x)
    g_n(i,2) = integral(@(z)g_int_arg(x(i),z), -Inf, Inf);
    G_n(i,2) = integral(@(z)G_int_arg(x(i),z), -Inf, Inf);
end
