function NEU=xyz84toNEU(xyz84,Pxyz84);%WGS84坐标系向站心直角坐标系转化
    lbh=xyz84tolbh(Pxyz84);
    l=lbh(1,1);
    b=lbh(1,2);
    CT=[-cos(l)*sin(b) -sin(l)*sin(b) cos(b);-sin(l) cos(l) 0;cos(l)*cos(b) sin(l)*cos(b) sin(b)];
    drtxyz84=[xyz84(:,1)-Pxyz84(1,1) xyz84(:,2)-Pxyz84(1,2) xyz84(:,3)-Pxyz84(1,3)];
    NEU=(CT*(drtxyz84'))';
end