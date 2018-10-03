%Statistical analysis for maximum sum rate for a collection of users in a fading channel with MRT beamforming.

rv_len = 10000;
%NOTE: l cannot be = 1!!!
% the sums used below will reduce the random vectors to scalars if l = 1
% (ie. the rows will be summed instead of the columns)
l = 2;  %group size
try
    assert(l>1);
catch me
    error(['Group size l (%i) out of range (l must be > 1). Error details: %s'],...
        l, me.message)
end
N = 4;          %num TX ants
tx_power = 1;   %total tx power in Watts
n = 6;          %size of candidate set

sigma_n = 0.1;  %noise variance


%cell array that contains a circularly symmetric normal rv that represents
% the channel for a given user and given antenna. each user has a vector
% of N RVs to describe the MIMO channel (one for each antenna)
h_arr = cell(n,1);
for arr_idx = 1:n
    h_arr{arr_idx,1} = (1/sqrt(2*N))*(randn(rv_len,l,1) + i*randn(rv_len,l,1));
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
            
    end
end

%model the sum rate as a normal distribution by central limit theorem.
% we use the properties of dependent random variables to come up with the
% pdf for the maximum sum rate accross all the sus groups.

%build up to evaluate the multi-var cdf, marginal pdf on
mu = mean_vec(1);
sigma = sqrt(cov_mat(1,1));
x = mu+(-2*sigma:sigma/10:8*sigma); %bounds of grid are estimated from first rv
x_grid = [];
for x_dim_idx = 1:num_groups
    x_grid = [x_grid x'];
end

cdf_mat = [[]];
pdf_mat = [[]];
for minus1_idx = 1:num_groups
    %keep grid size consistent
    x_grid_minus1 = x_grid; 
    x_grid_minus1(:,minus1_idx) = [];

    %conditional covariance
    cov_mat_minus1 = cov_mat;
    cov_mat_minus1(:,minus1_idx) = [];
    cov_mat_minus1(minus1_idx,:) = [];
    
    var_minus1 = cov_mat(minus1_idx,minus1_idx);
   
    cov_col_vec_minus1 = cov_mat(:,minus1_idx);
    cov_col_vec_minus1(minus1_idx) = [];
    
    cov_row_vec_minus1 = cov_mat(minus1_idx,:);
    cov_row_vec_minus1(minus1_idx) = [];

    cov_cdl_minus1 = cov_mat_minus1 - cov_col_vec_minus1*(cov_col_vec_minus1')/var_minus1;
   
    %conditional mean 
    mean_vec_minus1 = [mean_vec(1:minus1_idx-1)' mean_vec(minus1_idx+1:end)']';
    mean_minus1 = mean_vec(minus1_idx);
    for x_idx = 1:length(x_grid_minus1)
        mean_cdl_minus1 = mean_vec_minus1 + (x_grid(x_idx,1)-mean_minus1)*cov_col_vec_minus1/var_minus1;
        %multi-variant normal cdf with one element removed
        cdf_mat(x_idx,minus1_idx) = mvncdf(x_grid_minus1(x_idx,:),mean_cdl_minus1',cov_cdl_minus1);
    end
    %normal marginal pdf corresponding to the distribution that was removed from the multi-var cdf
    pdf_mat(:,minus1_idx) = normpdf(x_grid_minus1(:,1),mean_minus1,var_minus1); 
end

%according to R.B Arellano-Valle (2008) corollary 4, the pdf of the maximum is the sum of the product
% of the multivariant cdf and the marginal pdf. we expect a skew pdf.
max_pdf = sum((cdf_mat.*pdf_mat)')';
