sigma_n_sq = 0.1;   %Noise power
Ptx = 1;            %Tx power
N = 4;              %Number of transmit antennas
l = 4;              %Number of users in a group
n = 300;            %Number of available candidate users
g_max = 5;          %Max SINR realization value
SINR_rv_len = 1000; %Number of realizations in SINR RV
num_trials = 100000;%Number of trials used to estimate integral for MC integration
num_groups = nchoosek(n,l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section of code calculates the PDF for the SINR according to Feng et al (2014).
%   This calculation requires an integration on the dummy variable, nu. The result
%   of this integration is a PDF of the SINR that takes arguments gamma or g in this script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constant that is independent of nu
C_gamma = @(g)((sigma_n_sq.^N).*exp(-1.*N./g).*(N.^(l+N-1)))./((Ptx.^N).*gamma(N).*gamma(l-1).*(g.^2));
f_gamma_int_arg = @(g,v) (((1./g)-v).^(l-2)).*exp(N.*v-(N.*sigma_n_sq./(Ptx.*v)))./(v.^(N+1));
mu_int_arg = @(g) f_gamma(g);
f_gamma_est = zeros(1,SINR_rv_len);
g_vec = linspace(0.1,g_max,SINR_rv_len);
for est_idx = 1:length(f_gamma_est)
    g = g_vec(est_idx);
    int_res = integral(@(v)f_gamma_int_arg(g,v),0 ,1./g);
    f_gamma_est(est_idx) = C_gamma(g)*int_res;
end
%Check to see that the PDF sums to ~1
sum_pdf = trapz(g_vec,f_gamma_est);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section of code calculates the mean sum rate based on the SINR PDF previously
%   calculated. This mean value is calculated via trapezoid method, and via Monte
%   Carlo integration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Capacity terms and calculated mean value according to trapezoid method
log_term = log2(1+g_vec);
log_pdf_prod = log_term.*f_gamma_est;
mu_trap = l*trapz(g_vec,log_pdf_prod);

%Capacity terms and calculated variance value according to trapezoid method
log_pdf_prod = (log_term.^2).*f_gamma_est;
sigma_sq_trap = l*trapz(g_vec,log_pdf_prod);


%Monte Carlo integration.
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

mu = (l/num_trials)*sum(sum_arg);

f_cap = @(g) log2(1+g).^2;
sum_arg = f_cap(trunc_g_vals).*trunc_p_vals;

sigma_sq = (l/num_trials)*sum(sum_arg);
sigma = sqrt(sigma_sq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section of code calculates the PDF/CDF of the max amongst the various 
%   groups assuming the mean and variance values previously calculated. The
%   PDFs are assumed to be Gaussian, which holds for sufficiently large l by 
%   central limit theorem. The correlation value is also set to a value that
%   assumes that nchoosek(n,l) is sufficiently large. The method assumed here
%   comes from Tong's book The Multi-variate Normal Distribution (1990)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correlation Coefficient between groups:
%   Note that this assignment only holds for large enough binom(n,l)
rho = 1/l;


%Max order statistic for the un-correlated case
f_n = @(z) num_groups*(normcdf(z,0,1).^(num_groups-1)).*normpdf(z,0,1);
F_n = @(z) normcdf(z,0,1).^num_groups;

%Scale the normal distributions by the appropriate mean, variance
g_int_arg = @(x,z) (1/(sigma*sqrt(1-rho)))*f_n((((x-mu)/sigma) + sqrt(rho)*z)/(sqrt(1-rho))).*normpdf(z,0,1);
G_int_arg = @(x,z) (F_n((((x-mu)/sigma) + sqrt(rho)*z)/(sqrt(1-rho)))).*normpdf(z,0,1);


clear max_pdf_vec;
clear max_cdf_vec;
x = mu+(-8*sigma:sigma/30:8*sigma); %bounds could be variable here to whatever
max_pdf_vec(:,1) = x;
max_cdf_vec(:,1) = x;
for i = 1:length(x)
    max_pdf_vec(i,2) = integral(@(z)g_int_arg(x(i),z), -Inf, Inf);
    max_cdf_vec(i,2) = integral(@(z)G_int_arg(x(i),z), -Inf, Inf);
end
