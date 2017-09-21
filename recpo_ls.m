function [pos, El, sat_position,GDOP, basic_obs] = recpo_ls(obs,sats,time,Eph)
% RECPO_LS Computation of receiver position from pseudoranges
%          using ordinary least-squares principle

%Kai Borre 31-10-2001
%Copyright (c) by Kai Borre
%$Revision: 1.1 $  $Date: 2002/07/10  $
%     pos1 = zeros(4,1);
%     pos1(1:3) = [-2170110.089  4385047.662  4078168.881];
%     ELS = get_el(pos1,obs,sats,time,Eph);
%     idx=find(ELS<10);
%     obs(idx,:)=[];
%     sats(idx,:)=[];
v_light = 299792458;
dtr = pi/180;
m = size(obs,1);  % number of svs
el = zeros(m,1);
NOBD=length(find(sats(:)>40));
% identify ephemerides columns in Eph

for t = 1:m
    col_Eph(t) = find_eph(Eph,sats(t),time);
end

no_iterations = 6; %迭代次数
ps_corr = [];
sat_pos = [];
sat_position=zeros(m,1);
if NOBD>0&&m>NOBD  %判断是否为双系统
pos = zeros(5,1);
% pos(1:3) = [0 0 0];  %测站坐标赋予初值
pos(1:3) = [-2404075.5526   4510736.0036   3802326.9417];  %dorm测站坐标赋予初值
% pos(1:3) = [-2170122.144 	4385062.372	4078164.873];  %测站坐标赋予初值
for iter = 1:no_iterations
    A = [];
    omc = []; % observed minus computed observation
    ELS=[];%高度角
    for i = 1:m        
        k = col_Eph(i);
        if Eph(1,k)<=40
            tx_RAW = time - obs(i)/v_light;
            t0c = Eph(21,k);
            dt = check_t(tx_RAW-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_GPS = tx_RAW-tcorr;
            dt = check_t(tx_GPS-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_GPS = tx_RAW-tcorr;
            X = satpos(tx_GPS, Eph(:,k));
        else
            tx_RAW = time - obs(i)/v_light-14;%BD时间比GPS时间慢14s
            t0c = Eph(21,k);
            dt = check_t(tx_RAW-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_BD = tx_RAW-tcorr;
            dt = check_t(tx_BD-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_BD = tx_RAW-tcorr;
            X = satpos(tx_BD, Eph(:,k));
        end
        ELS=[ELS;solve_el(pos(1:3),X)];
   %X = satt(:,i);
   %for i=1:m
        
        if iter == 1
            traveltime = 0.072;
            Rot_X = X;
            trop = 0;
        else
            rho2 = (X(1)-pos(1))^2+(X(2)-pos(2))^2+(X(3)-pos(3))^2;
            traveltime = sqrt(rho2)/v_light;
            Rot_X = e_r_corr(traveltime,X);
            rho2 = (Rot_X(1)-pos(1))^2+(Rot_X(2)-pos(2))^2+(Rot_X(3)-pos(3))^2;          
            [az,el,dist] = topocent(pos(1:3,:),Rot_X-pos(1:3,:));                                                            
            if iter == no_iterations, El(i) = el; end
            trop = tropo(sin(el*dtr),0.0,1013.0,293.0,50.0,...  %这里用的是74年的一个对流层模型
                0.0,0.0,0.0);   
        end
        % subtraction of pos(4) corrects for receiver clock offset and
        % v_light*tcorr is the satellite clock offset
        %pos(5) 是相对于BD系统的接收机钟差
        if iter == no_iterations
            ps_corr = [ps_corr; obs(i)+v_light*tcorr-trop];
            sat_pos = [sat_pos; X'];
        end
        if i<=m-NOBD
        sat_position(i)=norm(Rot_X-pos(1:3))-v_light*tcorr;
        omc = [omc; obs(i)-norm(Rot_X-pos(1:3),'fro')-pos(4)+v_light*tcorr-trop];
        A = [A; (-(Rot_X(1)-pos(1)))/norm(Rot_X-pos(1:3))...
                (-(Rot_X(2)-pos(2)))/norm(Rot_X-pos(1:3))...
                (-(Rot_X(3)-pos(3)))/norm(Rot_X-pos(1:3)) 1 0];
        else
            sat_position(i)=norm(Rot_X-pos(1:3))-v_light*tcorr;
            omc = [omc; obs(i)-norm(Rot_X-pos(1:3),'fro')-pos(5)+v_light*tcorr-trop];
            A = [A; (-(Rot_X(1)-pos(1)))/norm(Rot_X-pos(1:3))...
                (-(Rot_X(2)-pos(2)))/norm(Rot_X-pos(1:3))...
                (-(Rot_X(3)-pos(3)))/norm(Rot_X-pos(1:3)) 0 1];
        end
    end % i
    W=zeros(m,m);
    for i=1:m
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
%     x = inv(A'*A)*(A'*omc);
    pos = pos+x;
%     sat_position=sat_position-obs;
    if iter == no_iterations, GDOP = sqrt(trace(inv(A'*A))); 
        %% two lines that solve an exercise on computing tdop
        % invm = inv(A'*A);
        % tdop = sqrt(invm(4,4))
    end
end % iter
else  %判断是否为双系统
pos = zeros(4,1);
% pos(1:3) = [0 0 0];  %测站坐标赋予初值
pos(1:3) = [-2404075.5526   4510736.0036   3802326.9417];  %dorm测站坐标赋予初值
% pos(1:3) = [-2403323.0859   4510529.6309   3803004.0969];  %torch测站坐标赋予初值
for iter = 1:no_iterations
    A = [];
    omc = []; % observed minus computed observation
    ELS=[];%高度角
    for i = 1:m        
        k = col_Eph(i);
        if Eph(1,k)<=40
            tx_RAW = time - obs(i)/v_light;
            t0c = Eph(21,k);
            dt = check_t(tx_RAW-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_GPS = tx_RAW-tcorr;
            dt = check_t(tx_GPS-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_GPS = tx_RAW-tcorr;
            X = satpos(tx_GPS, Eph(:,k));
        else
            tx_RAW = time - obs(i)/v_light-14;%BD时间比GPS时间慢14s
            t0c = Eph(21,k);
            dt = check_t(tx_RAW-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_BD = tx_RAW-tcorr;
            dt = check_t(tx_BD-t0c);
            tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
            tx_BD = tx_RAW-tcorr;
            X = satpos(tx_BD, Eph(:,k));
        end
        ELS=[ELS;solve_el(pos(1:3),X)];
   %X = satt(:,i);
   %for i=1:m
        
        if iter == 1
            traveltime = 0.072;
            Rot_X = X;
            trop = 0;
        else
            rho2 = (X(1)-pos(1))^2+(X(2)-pos(2))^2+(X(3)-pos(3))^2;
            traveltime = sqrt(rho2)/v_light;
            Rot_X = e_r_corr(traveltime,X);
            rho2 = (Rot_X(1)-pos(1))^2+(Rot_X(2)-pos(2))^2+(Rot_X(3)-pos(3))^2;          
            [az,el,dist] = topocent(pos(1:3,:),Rot_X-pos(1:3,:));                                                            
            if iter == no_iterations, El(i) = el; end
            trop = tropo(sin(el*dtr),0.0,1013.0,293.0,50.0,...
                0.0,0.0,0.0);   
        end
        % subtraction of pos(4) corrects for receiver clock offset and
        % v_light*tcorr is the satellite clock offset
        if iter == no_iterations
            ps_corr = [ps_corr; obs(i)+v_light*tcorr-trop];
            sat_pos = [sat_pos; X'];
        end
        sat_position(i)=norm(Rot_X-pos(1:3))-v_light*tcorr;
        omc = [omc; obs(i)-norm(Rot_X-pos(1:3),'fro')-pos(4)+v_light*tcorr-trop];
        A = [A; (-(Rot_X(1)-pos(1)))/norm(Rot_X-pos(1:3))...
                (-(Rot_X(2)-pos(2)))/norm(Rot_X-pos(1:3))...
                (-(Rot_X(3)-pos(3)))/norm(Rot_X-pos(1:3)) 1 ];
        end
    end % i
    W=zeros(m,m);
    for i=1:m
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
    pos = pos+x;
    sat_position=sat_position-obs;
    if iter == no_iterations, GDOP = sqrt(trace(inv(A'*A))); 
        %% two lines that solve an exercise on computing tdop
        % invm = inv(A'*A);
        % tdop = sqrt(invm(4,4))
    end
end % iter
basic_obs = [sat_pos ps_corr];

%%%%%%%%%%%%%%%%%%%%%  recpo_ls.m  %%%%%%%%%%%%%%%%%%%%%
