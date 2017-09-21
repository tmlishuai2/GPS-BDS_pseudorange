function [ eph ] = findeph( t, sv, nav )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
index = find(nav(:,1) == sv);
time_n = cal2gpst(nav(index,2:7));
time_n = time_n(:,2);
[m,i] = min(abs(time_n - t));
eph = nav(index(i),:);
end

