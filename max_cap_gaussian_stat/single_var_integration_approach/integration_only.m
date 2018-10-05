l = 4;
n = 100;

num_groups = nchoosek(n,l);

mu = 1.19;
rho = 1/l;


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
