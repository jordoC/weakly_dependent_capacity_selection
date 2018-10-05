%Proof of concept script to investigate results given by Tong (1990)

n = 10;     %number of normal distributions to sample from
sigma = 1;  %variance of each dist
rho = .99;  %correlation between distribtuions
mu = 2.5;   %mean of each dist

f_n = @(z) n*(normcdf(z,0,1).^(n-1)).*normpdf(z,0,1);
F_n = @(z) normcdf(z,0,1);

g_int_arg = @(x,z) (1/(sigma*sqrt(1-rho)))*f_n((((x-mu)/sigma) + sqrt(rho)*z)/(sqrt(1-rho))).*normpdf(z,0,1);
G_int_arg = @(x,z) (F_n((((x-mu)/sigma) + sqrt(rho)*z)/(sqrt(1-rho))).^n).*normpdf(z,0,1);


clear g_n;
clear G_n;
x = 0:0.1:10;
g_n(:,1) = x;
G_n(:,1) = x;
for i = 1:length(x)
    g_n(i,2) = integral(@(z)g_int_arg(x(i),z), -Inf, Inf);
    G_n(i,2) = integral(@(z)G_int_arg(x(i),z), -Inf, Inf);
end
