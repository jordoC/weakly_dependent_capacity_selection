function [integral_res] = mc_integral(f,a,b,err)
    %integrate function (f) using monte carlo method
    N = 1/(err^2);
    x_rnd = (b-a)*rand(1,N)+a;
    f_rnd = f(x_rnd);
    integral_res = ((b-a)/numel(x_rnd))*sum(f_rnd);
end
%%function [integral_dx] = mc_integral( f,a,b )
%f = @(x) x;
%a = 1;
%b = 2;
%    %integrate function (f) using monte carlo method
%    x2=linspace(a,b,1000);
%    syms z % zero vector holder to find max y value
%    z = zeros(size(x2));
%    z = f(x2);
%    y = f(b).*rand(1,1000);
%    x = rand(1,1000);
%    h=0; % counters
%    n= 0;
%    %if you want to see visual representation just un-commnet plot lines
%    clf;
%    plot(x2,z); hold on;
%    plot(x,y,'x')
%    count = 0;
%    for k=1:numel(x);
%        if y(k) <= exp(x(k))+1;
%            count= count +1;
%        end
%    end
%    count
%    integral_dx = count/numel(x) * max(z) * (b-a)
%%end
