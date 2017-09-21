function satp = satpos(t,eph);
%SATPOS   Calculation of X,Y,Z coordinates at time t
%         for given ephemeris eph

%Kai Borre 04-09-96
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 1997/09/26  $

%eph为单条星历数据
GM_GPS = 3.986005e14;           %GPS地球引力常数
GM_BD = 3.986004418e14;         %BD地球引力常数
% parameter m^3/s^2
Omegae_dot_GPS = 7.2921151467e-5; %GPS地球旋转速率
Omegae_dot_BD = 7.292115e-5;  %BD地球旋转速率
pi=3.1415926535898;  %pi取值，GPS与BD一样
%  Units are either seconds, meters, or radians
%  Assigning the local variables to eph
svprn   =   eph(1);
af2     =   eph(10);
M0      =   eph(14);
roota   =   eph(18);
deltan  =   eph(13);
ecc     =   eph(16);
omega   =   eph(25);
cuc     =   eph(15);
cus     =   eph(17);
crc     =  eph(24);
crs     =  eph(12);
i0      =  eph(23);
idot    =  eph(27);
cic     =  eph(20);
cis     =  eph(22);
Omega0  =  eph(21);
Omegadot=  eph(26);
toe     =  eph(19);
af0     =  eph(8);
af1     =  eph(9);
toc     =  eph(2:7);
toc = cal2gpst(toc);
toc = toc(:,2);
% Procedure for coordinate calculation
% t = 565760.87291000085
%计算GPS卫星位置
if floor(svprn/100)==71
    A = roota*roota;
    tk = check_t(t-toe);
    n0 = sqrt(GM_GPS/A^3);
    n = n0+deltan;
    M = M0+n*tk;
    M = rem(M+2*pi,2*pi);
    E = M;
    for i = 1:10
        E_old = E;
        E = M+ecc*sin(E);
        dE = rem(E-E_old,2*pi);
        if abs(dE) < 1.e-16
            break;
        end
    end
    E = rem(E+2*pi,2*pi);
    v = atan2(sqrt(1-ecc^2)*sin(E), cos(E)-ecc);
    phi = v+omega;
    phi = rem(phi,2*pi);
    u = phi              + cuc*cos(2*phi)+cus*sin(2*phi);
    r = A*(1-ecc*cos(E)) + crc*cos(2*phi)+crs*sin(2*phi);
    i = i0+idot*tk       + cic*cos(2*phi)+cis*sin(2*phi);
    Omega = Omega0+(Omegadot-Omegae_dot_GPS)*tk-Omegae_dot_GPS*toe;
    Omega = rem(Omega+2*pi,2*pi);
    x1 = cos(u)*r;
    y1 = sin(u)*r;
    satp(1,1) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
    satp(2,1) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
    satp(3,1) = y1*sin(i);
elseif mod(svprn,100)>=1 && mod(svprn,100)<=5
        A = roota*roota;
        tk = check_t(t-toe);
        n0 = sqrt(GM_BD/A^3);
        n = n0+deltan;
        M = M0+n*tk;
        M = rem(M+2*pi,2*pi);
        E = M;
        for i = 1:10
            E_old = E;
            E = M+ecc*sin(E);
            dE = rem(E-E_old,2*pi);
            if abs(dE) < 1.e-16
                break;
            end
        end
        E = rem(E+2*pi,2*pi);
        v = atan2(sqrt(1-ecc^2)*sin(E), cos(E)-ecc);
        phi = v+omega;
        phi = rem(phi,2*pi);
        u = phi              + cuc*cos(2*phi)+cus*sin(2*phi);
        r = A*(1-ecc*cos(E)) + crc*cos(2*phi)+crs*sin(2*phi);
        i = i0+idot*tk       + cic*cos(2*phi)+cis*sin(2*phi);
        Omega = Omega0+Omegadot*tk-Omegae_dot_BD*toe;
        Omega = rem(Omega+2*pi,2*pi);
        x1 = cos(u)*r;
        y1 = sin(u)*r;
        satp1(1,1) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
        satp1(2,1) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
        satp1(3,1) = y1*sin(i);
        Rx = [1 0 0
              0 cos(-5*pi/180) sin(-5*pi/180)
              0 -sin(-5*pi/180) cos(-5*pi/180)];
        Rz = [cos(Omegae_dot_BD*tk) sin(Omegae_dot_BD*tk) 0
              -sin(Omegae_dot_BD*tk) cos(Omegae_dot_BD*tk) 0
              0 0 1]; 
        temp = Rx*satp1;
        satp = Rz*temp;
else
            A = roota*roota;
            tk = check_t(t-toe);
            n0 = sqrt(GM_BD/A^3);
            n = n0+deltan;
            M = M0+n*tk;
            M = rem(M+2*pi,2*pi);
            E = M;
            for i = 1:10
                E_old = E;
                E = M+ecc*sin(E);
                dE = rem(E-E_old,2*pi);
                if abs(dE) < 1.e-16
                    break;
                end
            end
            E = rem(E+2*pi,2*pi);
            v = atan2(sqrt(1-ecc^2)*sin(E), cos(E)-ecc);
            phi = v+omega;
            phi = rem(phi,2*pi);
            u = phi              + cuc*cos(2*phi)+cus*sin(2*phi);
            r = A*(1-ecc*cos(E)) + crc*cos(2*phi)+crs*sin(2*phi);
            i = i0+idot*tk       + cic*cos(2*phi)+cis*sin(2*phi);
            Omega = Omega0+(Omegadot-Omegae_dot_BD)*tk-Omegae_dot_BD*toe;
            Omega = rem(Omega+2*pi,2*pi);
            x1 = cos(u)*r;
            y1 = sin(u)*r;
            satp(1,1) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
            satp(2,1) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
            satp(3,1) = y1*sin(i);   
end
end
%%%%%%%%% end satpos.m %%%%%%%%%
