function pn = premium_notional(n, distribution)  
    values=(0:(size(distribution,1)-1))';
    pn=n-sum(min(n,values).*distribution);
end