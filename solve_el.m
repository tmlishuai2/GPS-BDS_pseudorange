function EL=solve_el(pos,sats)%3*1æÿ’Û
    p0_lbh=xyz84tolbh(pos');
    CT=[-cos(p0_lbh(1))*sin(p0_lbh(2)) -sin(p0_lbh(1))*sin(p0_lbh(2)) cos(p0_lbh(2));
        -sin(p0_lbh(1)) cos(p0_lbh(1)) 0;
        cos(p0_lbh(1))*cos(p0_lbh(2)) sin(p0_lbh(1))*cos(p0_lbh(2)) sin(p0_lbh(2))];
    drtxyz84=sats-pos;
    NEU=(CT*(drtxyz84));
    R=sqrt(NEU(1)^2+NEU(2)^2+NEU(3)^2);
    EL=asin(NEU(3)/R);%ª°∂»
end