function std_1=std_1(v);
    n=length(v);
    m=0;
    for i=1:n
        m=m+v(i)^2;
    end;
    std_1=sqrt(m/n);
end