function [El,az,XXX] = get_elaz(pos,obs,sats,time,Eph)     
%输出此历元各个卫星相对于测站的高度角和方位角以及各个卫星的位置
v_light = 299792458;
dtr = pi/180;
m = size(obs,1);  % number of svs
% identify ephemerides columns in Eph
for t = 1:m
    col_Eph(t) = find_eph(Eph,sats(t),time);
end

% preliminary guess for receiver position and receiver clock offset
for i = 1:m 
    k = col_Eph(i);
    tx_RAW = time - obs(i)/v_light;
    %***********调整GPST与BDT***************
    if  Eph(1,k) > 40
        tx_RAW = tx_RAW - 14;
    end
    %***************************************
    t0c = Eph(21,k);
    dt = check_t(tx_RAW-t0c);
    tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
    tx_GPS = tx_RAW-tcorr;
    dt = check_t(tx_GPS-t0c);
    tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
    tx_GPS = tx_RAW-tcorr;
    X = satpos(tx_GPS, Eph(:,k));
    XXX(:,i)=X;
    [az(i),El(i),dist] = topocent(pos(1:3,:),X-pos(1:3,:));                                                              
end % iter

%%%%%%%%%%%%%%%%%%%%%  get_el.m  %%%%%%%%%%%%%%%%%%%%%

