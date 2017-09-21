clear all;
clc;
global f0 f1 f2 f5 v_light b1 b2 b3;
f0 = 10.23*10^6;
f1 = 154*f0;
f2 = 120*f0;
f5 = 115*f0;
v_light = 299792458;
b1 = 1561.098*10^6;
b2 = 1207.14*10^6;
b3 = 1268.52*10^6;
% [nhdr, nav] = read_rinex_nav('brdm1130.16p');
% [ohdr, obs] = read_rinex_obs('B0051130.16O');
%obs文件读取速度比较慢，因此这里直接选择导入已经读取好的矩阵文件
load('b0051130.16o.mat');
load('brdm1130.16p.mat');

NoObs_types1 = ohdr.nobs;%观测类型数
%计算有多少个历元，每个历元多少个观测值
gpst_o = cal2gpst(obs(:,1:6));
gpst_o = gpst_o(:,2);
nepch_o = 1;
nsat_o = [];
for k = 2:size(gpst_o)
    if gpst_o(k-1) < gpst_o(k)
        nepch_o = nepch_o + 1;
        nsat_o = [nsat_o; k-1];
    end
end
nsat_o1 = [nsat_o(2:end); size(gpst_o,1)];
nsat_o2 = nsat_o1 - nsat_o; 
nsat_o = [nsat_o(1); nsat_o2];
NoC1C = 9;
NoC2P = 13;
NoC5I = 17;
NoC2I = 9;
NoC7I = 13;
NoC6I = 17;

%将nav中的GPS与BDS的信息提取出来
index = find((floor(nav(:,1)/100)==71)+(floor(nav(:,1)/100)==67));
nav = nav(index,:);
gpst_n = cal2gpst(nav(:,2:7));
gpst_n = gpst_n(:,2);

Sta1 = [-1351114.3841   5603358.2193   2726446.5904]';    %测站的精确坐标，由rtklib长时间解算得到，解算中并未用到，用于比较结果
        
Pos = [];
% sat_El = zeros(10,qend/interval);
GDOP=[];

epoch = 1:420;   %epoch是用户给定的需要求解的历元时刻
count = size(epoch,2); 

start = 1;
%每个历元求解一次
for q = epoch
    dtr = pi/180;
    m = nsat_o(q);  % number of svs
    

    El = [];
    ps_corr = [];
    sat_pos = [];
    obsq = obs(start:start+nsat_o(q)-1,:);
    time = gpst_o(start);
    sv = obsq(:,8);
    NOBD=length(find(floor(sv/100)==67));
    pos = zeros(5,1);
    pos(1:3) = [-1341121.07800000,5613385.02650000,2736460.75330000];
    no_iterations = 6;  
    
    for iter = 1:no_iterations
        A = [];
        omc = []; % observed minus computed observation
        ELS=[];
        sat_position = [];
%        el = zeros(m,1);
        %找time和sv对应的星历数据
        for i = 1:size(sv,1)
            if floor(sv(i)/100) == 71
                obsqi = ( f2^2*obsq(i,13) - f1^2*obsq(i,9) )/( f2^2 - f1^2 );      %w无电离层组合
                if isnan(obsqi)
                    continue;
                end
                eph = findeph(time,sv(i),nav);
                tx_RAW = time - obsqi/v_light;  %
                t0c = cal2gpst(eph(2:7));
                t0c = t0c(2);
                dt = check_t(tx_RAW-t0c);
                tcorr = (eph(10)*dt + eph(9))*dt + eph(8);
                tx_GPS = tx_RAW-tcorr;
                dt = check_t(tx_GPS-t0c);
                tcorr = (eph(10)*dt + eph(9))*dt + eph(8);
                tx_GPS = tx_RAW-tcorr;
                X = satpos(tx_GPS, eph);
            else
                obsqi = ( b2^2*obsq(i,13) - b1^2*obsq(i,9) )/( b2^2 - b1^2 );      %w无电离层组合
                if isnan(obsqi)
                    continue;
                end
                eph = findeph(time,sv(i),nav);
                tx_RAW = time - obsqi/v_light -14; %BDS比GPS慢14s
                t0c = cal2gpst(eph(2:7));
                t0c = t0c(2);
                dt = check_t(tx_RAW-t0c);
                tcorr = (eph(10)*dt + eph(9))*dt + eph(8);
                tx_BDS = tx_RAW-tcorr;
                dt = check_t(tx_BDS-t0c);
                tcorr = (eph(10)*dt + eph(9))*dt + eph(8);
                tx_BDS = tx_RAW-tcorr;
                X = satpos(tx_BDS, eph);
            end
            ELS=[ELS;solve_el(pos(1:3),X)];
       %X = satt(:,i);
       %for i=1:m

            if iter == 1
                traveltime = 0.072;
                Rot_X = X;
                trop = 0;
            else
                traveltime = norm(X-pos(1:3))/v_light;
                Rot_X = e_r_corr(traveltime,X);
%                 rho2 = (Rot_X(1)-pos(1))^2+(Rot_X(2)-pos(2))^2+(Rot_X(3)-pos(3))^2;          
                [az,el,dist] = topocent(pos(1:3,:),Rot_X-pos(1:3,:));                                                            
                if iter == no_iterations
                    El = [El; el]; 
                end
                trop = tropo(sin(el*dtr),0.0,1013.0,293.0,50.0,...  %这里用的是74年的一个对流层模型
                    0.0,0.0,0.0);   
            end
            % subtraction of pos(4) corrects for receiver clock offset and
            % v_light*tcorr is the satellite clock offset
            %pos(5) 是相对于BD系统的接收机钟差
            if iter == no_iterations
                ps_corr = [ps_corr; obsqi+v_light*tcorr-trop];
                sat_pos = [sat_pos; X'];
            end
            if floor(sv(i)/100) == 71
                sat_position(i)=norm(Rot_X-pos(1:3))-v_light*tcorr;
                omc = [omc; obsqi-norm(Rot_X-pos(1:3),'fro')-pos(4)+v_light*tcorr-trop];
                A = [A; (-(Rot_X(1)-pos(1)))/norm(Rot_X-pos(1:3))...
                        (-(Rot_X(2)-pos(2)))/norm(Rot_X-pos(1:3))...
                        (-(Rot_X(3)-pos(3)))/norm(Rot_X-pos(1:3)) 1 0];
            else
                sat_position(i)=norm(Rot_X-pos(1:3))-v_light*tcorr;
                omc = [omc; obsqi-norm(Rot_X-pos(1:3),'fro')-pos(5)+v_light*tcorr-trop];
                A = [A; (-(Rot_X(1)-pos(1)))/norm(Rot_X-pos(1:3))...
                    (-(Rot_X(2)-pos(2)))/norm(Rot_X-pos(1:3))...
                    (-(Rot_X(3)-pos(3)))/norm(Rot_X-pos(1:3)) 0 1];
             end
        end % i
        W=eye(size(A,1));
        for i=1:size(A,1)
            W(i,i)=sin(ELS(i))^2;
        end

    %      for i=1:m-5
    %          W(i,i)=sin(ELS(i))^2*2;
    %      end
    %      for i=m-5+1:m
    %          W(i,i)=sin(ELS(i))^2;
    %      end

    %     x = inv(A'*A)*(A'*omc);
    
        x = inv(A'*W*A)*(A'*W*omc);
%         x = A\omc;
    %     x = inv(A'*A)*(A'*omc);
        pos = pos+x;
    %     sat_position=sat_position-obs;
        if iter == no_iterations
            gdop = sqrt(trace(inv(A'*A))); 
            % two lines that solve an exercise on computing tdop
            % invm = inv(A'*A);
            % tdop = sqrt(invm(4,4))
        end
    end % iter
    start = start + nsat_o(q);
    Pos = [Pos pos];
    GDOP = [GDOP gdop];
end   

me = mean(Pos,2);
agdop = mean(GDOP);
clockbias_g = mean(Pos(4,:));
clockbias_c = mean(Pos(5,:));
fprintf('\nMean GDOP: %12.3f\n', agdop)
fprintf('\nMean receiver clock bias in GPS: %12.3f\n', clockbias_g);
fprintf('\nMean receiver clock bias in BDS: %12.3f\n', clockbias_c);
fprintf('\nMean Position as Computed From %d Epochs:',count);
fprintf('\n\nAverage position：\nX: %12.3f  Y: %12.3f  Z: %12.3f\n', me(1,1), me(2,1), me(3,1));
fprintf('\n\nPositioning error:\ndeltaX: %12.3f  deltaY: %12.3f  deltaZ: %12.3f\n', me(1,1)-Sta1(1), me(2,1)-Sta1(2), me(3,1)-Sta1(3));   
        
%画位置误差散点图        
figure
plot((Pos(1:3,:)-Sta1*ones(1,count))','x')
%axis([0 420 -50 50])
legend('DeltaX','DeltaY','DeltaZ')
title('Change of WGS-84 coordinate positioning error','fontsize',12)
xlabel('Epoch（Interval:1s）','fontsize',12)
ylabel('Change of coordinate/m','fontsize',12)
set(gca,'fontsize',12)
  
        
%显示最值，均值和方差
MAX_xyz=zeros(3,4);
fprintf('\nStatistics of positioning in xyz(m)');
fprintf('\n     max       min      mean      std\n');
for i =1:3
    max_v = max(Pos(i,:)-Sta1(i,1)*ones(1,count));
    min_v = min(Pos(i,:)-Sta1(i,1)*ones(1,count));
    mean_v = mean(Pos(i,:)-Sta1(i,1)*ones(1,count));
    std_v = std_1(Pos(i,:)-Sta1(i,1)*ones(1,count));
    fprintf('\n%8.3f  %8.3f  %8.3f  %8.3f\n', max_v, min_v, mean_v, std_v);
    MAX_xyz(i,1)=max_v;
    MAX_xyz(i,2)=min_v;
    MAX_xyz(i,3)=mean_v;
    MAX_xyz(i,4)=std_v;
end

% 画NEU图
NEU=zeros(3,count);
for i=1:count
     NEU(:,i)=xyz84toNEU(Pos(1:3,i)',Sta1(1:3)')';
end
figure
plot(NEU','.','markersize',5)
%axis([0 420 -50 50])
legend('N','E','U')
title('Change of topocentric coordinate','fontsize',12)
xlabel('Epoch（Interval:1s）','fontsize',12)
ylabel('Change of coordinate/m','fontsize',12)
set(gca,'fontsize',12)


%显示NEU的最值，均值，方差 
MAX_neu=zeros(3,4);
fprintf('\nStatistics of positioning coordinate in NEU(m)\n   max       min      mean      std\n');
for i =1:3
    max_NEU = max(NEU(i,:));
    min_NEU = min(NEU(i,:));
    mean_NEU = mean(NEU(i,:));
    std_NEU = std_1(NEU(i,:));
    fprintf('\n%8.3f  %8.3f  %8.3f  %8.3f\n', max_NEU, min_NEU, mean_NEU, std_NEU);
    MAX_neu(i,1)=max_NEU;
    MAX_neu(i,2)=min_NEU;
    MAX_neu(i,3)=mean_NEU;
    MAX_neu(i,4)=std_NEU;
end

%画GDOP值变化
figure
plot(GDOP','.','markersize',5)
%axis([0 420 0 10])
legend('GDOP')
title('Change of GDOP','fontsize',12)
xlabel('Epoch（Interval:1s）','fontsize',12)
ylabel('GDOP','fontsize',12)
set(gca,'fontsize',12)

%画接收机钟差变化
figure
plot(Pos(4,:)/299792458*10^9,'.','markersize',5)
axis([0 420 (min(Pos(4,:))-10)/299792458*10^9 (max(Pos(4,:))+10)/299792458*10^9])
title('Change of receiver clock error in GPS','fontsize',12)
xlabel('Epoch（Interval:1s）','fontsize',12)
ylabel('Receiver clock error/ns','fontsize',12)
set(gca,'fontsize',12)

Sta1=me;
figure
plot((Pos(1:3,:)-Sta1(1:3,1)*ones(1,count))','.')
%axis([0 420 -50 50])
legend('DeltaX','DeltaY','DeltaZ')
title('Change of WGS-84 positioning coordinate','fontsize',12)
xlabel('Epoch（Interval:1s）','fontsize',12)
ylabel('Change of coordinate/m','fontsize',12)
set(gca,'fontsize',12)


        
        
        
        
        
        

