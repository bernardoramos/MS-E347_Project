function [Lt, Nt] = get_distribution(params, T, N, truncation)
    % truncation between 0 and 1
    u=zeros(N,2);  
    u(:,2)=2*i*pi*(0:N-1)'/N;
    
    [alpha, beta]=ab(u, params, T, 12);
    lambda0=params(4);
    transform=exp(alpha+beta*lambda0);
    distr=ifft(transform);  
    nt = floor(N*truncation);
    distr=distr(1:nt);
    distr=distr/size(distr,2);
    distr=abs(distr);
    Nt=distr;
    
    u=zeros(N,2);  
    u(:,1)=2*i*pi*(0:N-1)'/N;
    
    [alpha, beta]=ab(u, params, T, 12);
    lambda0=params(4);
    transform=exp(alpha+beta*lambda0);
    distr=ifft(transform);
    nt = floor(N*truncation);
    distr=distr(1:nt);
    distr=distr/size(distr,2);
    distr=abs(distr);
    Lt=distr;
    
end