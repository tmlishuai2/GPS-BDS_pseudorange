function dadi=xyz84tolbh(xyz);%WGS-84坐标系向大地坐标系转换
    x=xyz(:,1);
    y=xyz(:,2);
    z=xyz(:,3);
    m=length(x);
    a = 6378137;
    f = 1/298.257223563;
    e2 = 2*f-f^2;
    for i=1:m
        l(i,1)=atan(y(i,1)/x(i,1));
        if l(i,1)<0
            l(i,1)=l(i,1)+pi;
        end
        b(i,1)=atan((z(i,1)/sqrt(x(i,1)^2+y(i,1)^2))*(1+e2));
        b1(i,1)=0;
        while abs(b(i,1)-b1(i,1))>0.1^10
            b1(i,1)=b(i,1);
            n=a/sqrt(1-e2*(sin(b1(i,1)))^2);
            h1(i,1)=z(i,1)/sin(b1(i,1))-n*(1-e2);
            b(i,1)=atan(z(i,1)*(n+h1(i,1))/(sqrt((x(i,1)^2+y(i,1)^2))*(n*(1-e2)+h1(i,1))));
        end
    n=a/sqrt(1-e2*(sin(b(i,1)))^2);
    h(i,1)=z(i,1)/sin(b(i,1))-n*(1-e2);
    end
    dadi=[l b h];
end