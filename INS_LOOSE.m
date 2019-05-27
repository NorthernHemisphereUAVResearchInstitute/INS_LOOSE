%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         指北方位捷联式惯性导航系统 + GPS反馈松组合        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%导航坐标系为北东地坐标系
clear all;clc;
%定义常量
global S;
S.PI = 3.14159265358979323846 ;
S.We =0.00007292115;
S.Re =6378137.0;
S.e  =0.003352813;
S.Con_d2r    =1.745329251994330e-002;
S.Con_r2d    =5.729577951308232e+001;  %radian to degree transform coefficient /*/
S.Con_dh2rs  =4.848136811095360e-006;  %d/h to rad/s */

load coreTurn.mat;
GPSstop = fix(GPS(1,end));
inittime = 0.0;
%%
%初始化全局变量%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LI,I] = LIInit(inittime);
%%
while(LI.flagRun==true) 
%     [I.nRn,I.nRe,I.gn,I.Wie,I.Wep] = EarthCalc(LI.latitude,LI.height,LI.vn);
    [I.nRn,I.nRe,I.gn,I.Wie,I.Wep] = EarthCalc(40.006403*S.PI/180,50,LI.vn);
    if(MMQ(1,2*LI.indexINS)>=GPS(1,LI.indexGPS) && MMQ(1,2*LI.indexINS)<=GPSstop)
        [LI,I] = LIKalman(LI,I,GPS);        
        LI.indexGPS = LI.indexGPS +1;
        LI.indexGPS
    elseif(MMQ(1,2*LI.indexINS)>GPSstop)
       LI.flagRun = false;
    end
    %惯导解算
    [LI,I] = INSCycle2(LI,I,MMQ);
    %保存数据
    Save_pos(:,LI.indexINS)=[LI.latitude;LI.longitude;LI.height];
    Save_vn(:,LI.indexINS)=LI.vn;
    Save_phi(:,LI.indexINS)=DCM2Euler(I.Cnb,LI.gama);
    LI.sita=Save_phi(1,LI.indexINS);LI.gama=Save_phi(2,LI.indexINS);LI.psi=Save_phi(3,LI.indexINS);
    LI.indexINS = LI.indexINS + 1;
%     LI.indexINS
end
%%
save('all.mat');

% plot result
t=[1:length(Save_pos)]*2*I.ts;

figure(3),clf(figure(3));propedit(figure(3));set((figure(3)),'Color','w');                              
subplot(331), plot(t,(Save_pos(1,:)-Save_pos(1,1)),'b'); ylabel('latitude [m]'),grid on;                              
subplot(332), plot(t,(Save_pos(2,:)-Save_pos(2,1)),'b'); ylabel('longitude [m]'),grid on;                             
subplot(333), plot(t,(Save_pos(3,:)-Save_pos(3,1)),'b'); ylabel('height [m]'),grid on;
subplot(334), plot(t,Save_vn(1,:),'b'); ylabel('VN [m/s]'),grid on;                              
subplot(335), plot(t,Save_vn(2,:),'b'); ylabel('VE [m/s]'),grid on;                             
subplot(336), plot(t,Save_vn(3,:),'b'); ylabel('VD [m/s]'),grid on; 
subplot(337), plot(t,Save_phi(1,:)*180/S.PI,'b'); ylabel('pitch [deg]'),grid on;                              
subplot(338), plot(t,Save_phi(2,:)*180/S.PI,'b'); ylabel('roll [deg]'),grid on;xlabel('Time [sec]');                               
subplot(339), plot(t,Save_phi(3,:)*180/S.PI,'b'); ylabel('azimuth [deg]'),grid on; 

