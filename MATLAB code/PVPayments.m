function D = PVPayments(tranche, params, T, N, truncation, r)
    Ku=tranche(1);
    Kd=tranche(2);
    n=tranche(3);
       
    D=0;
    g=24;
    grille = linspace(0, T, g);
    for i=1:length(grille)
        s=grille(i);
        [d_tmp, ~ ] =get_distribution(params, s, N, truncation);
        values=.6*(0:length(d_tmp)-1)';
        e_tmp=sum((max(values - Kd * n, 0) - max(values - Ku * n, 0)).* d_tmp);
        D=D+r*exp(-r*s)*e_tmp*1/g;
    end
    D=D+exp(-r*T)*e_tmp;
end