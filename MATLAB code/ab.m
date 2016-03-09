function [a b] = ab(u, params, T, N)

% params = c, kappa, delta, lambda0
% N number of steps
% T maturity
% Initial condition : 0

h=T/N;
c=params(1);
kappa=params(2);
delta=params(3);

n=size(u,1);
tmp_b = zeros(n, 2*N+1);
h2=T/(2*N);

partial_b = @(y) (-kappa*y-1+exp(.6*(delta*y+u(:,1))+u(:,2)));

for i=1:(2*N)
    k1=partial_b(tmp_b(:,i));
    k2=partial_b(tmp_b(:,i)+h2/2*k1);
    k3=partial_b(tmp_b(:,i)+h2/2*k2);
    k4=partial_b(tmp_b(:,i)+h2*k3);
    
    tmp_b(:,i+1)=tmp_b(:,i)+h2/6*(k1+2*k2+2*k3+k4);
end

b=tmp_b(:,end);

tmp_a = zeros(n, N+1);

for i=1:N
    k1=kappa*c*tmp_b(:,2*i-1);
    k2=kappa*c*tmp_b(:,2*i);
    k3=kappa*c*tmp_b(:,2*i);
    k4=kappa*c*tmp_b(:,2*i+1);
    
    tmp_a(:,i+1)=tmp_a(:,i)+h/6*(k1+2*k2+2*k3+k4);
end

a=tmp_a(:,end);

end