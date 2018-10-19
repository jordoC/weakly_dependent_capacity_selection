sigma_n_sq = 0.1;
Ptx = 1;
N = 4;
l = 4;


C_gamma = @(g)((sigma_n_sq.^N).*exp(-1.*N./g).*(N.^(l+N-1)))./((Ptx.^N).*gamma(N).*gamma(l-1).*(g.^2));
f_gamma_int_arg = @(g,v) (((1./g)-v).^(l-2)).*exp(N.*v-(N.*sigma_n_sq./(Ptx.*v)))./(v.^(N+1));
%f_gamma = @(g) C_gamma(g)*integral(@(v)f_gamma_int_arg(g,v),0 ,1./g);
%f_gamma(1)
%mu_int_arg = @(g) log2(1+g).*f_gamma(g);
mu_int_arg = @(g) f_gamma(g);
%mu  = integral(@(g)mu_int_arg(g),0,Inf)
%f_n = @(z) num_groups*(normcdf(z,0,1).^(num_groups-1)).*normpdf(z,0,1);

%g_rnd = g_ub*rand(1,num_trials);
%v_rnd = rand(1,num_trials);
%v_rnd = v_rnd./g_rnd;
%
%
%f_rnd = f_gamma_int_arg(g_rnd,v_rnd);
%f_rnd_pdf = f_rnd./g_rnd;
%integral_res = (1/num_trials)*sum(f_rnd_pdf);

f_gamma_est = zeros(1,1000);
g_max =5; 
g_vec = linspace(0.1,g_max,1000);
num_trials = 100000;
for est_idx = 1:length(f_gamma_est)
    g = g_vec(est_idx);
    int_res = integral(@(v)f_gamma_int_arg(g,v),0 ,1./g);
    f_gamma_est(est_idx) = C_gamma(g)*int_res;
end


sum_pdf = trapz(g_vec,f_gamma_est);

log_term = log2(1+g_vec);
log_pdf_prod = log_term.*f_gamma_est;
mu_trap = l*trapz(g_vec,log_pdf_prod)

log_pdf_prod = (log_term.^2).*f_gamma_est;
sigma_trap = l*trapz(g_vec,log_pdf_prod)


pdf_gamma_est = [g_vec' f_gamma_est'];

rnd_g_vals = g_max*rand(1,num_trials);
trunc_g_vals = zeros(1,num_trials);
trunc_p_vals = zeros(1,num_trials);

for rnd_g_idx = 1:length(rnd_g_vals)
    rnd_g_val = rnd_g_vals(rnd_g_idx);
    [g_diff g_idx] = min(abs(pdf_gamma_est(:,1) - rnd_g_val));
    trunc_g_vals(rnd_g_idx) = pdf_gamma_est(g_idx,1);
    trunc_p_vals(rnd_g_idx) = pdf_gamma_est(g_idx,2);
end

f_cap = @(g) log2(1+g);
sum_arg = f_cap(trunc_g_vals).*trunc_p_vals;

mu = l*(g_max/num_trials)*sum(sum_arg)

f_cap = @(g) log2(1+g).^2;
sum_arg = f_cap(trunc_g_vals).*trunc_p_vals;

sigma = l*(g_max/num_trials)*sum(sum_arg)

