function S=model(tranche, params, T, N, truncation, r)
    Ku=tranche(1);
    Kd=tranche(2);
    n=tranche(3);
    is_upfront=tranche(4);
    
    [~, Nt] = get_distribution(params, T, N, truncation);
    D = PVPayments(tranche, params, T, N, truncation, r);
    
    if (is_upfront)
        K=n*(Ku-Kd);
        S=D/K;
    else
        sum=0;
        premium_dates=0:.25:(T-.25);
        cm=.25;
        for i=1:length(premium_dates)
            tm=premium_dates(i);
            sum=sum+exp(-r*tm)*cm*premium_notional(tm, Nt);
        end
        S=D/sum;
    end
end