clc;
clear all;

DATA=load('C:\Users\LuoHongfu\Desktop\INS_LOOSE\pull_test_13.csv');%�����ļ�����

first_start=1;
first_end=length(DATA(:,1));

distance=DATA(first_start:first_end,1)-DATA(first_start,1);
time=DATA(first_start:first_end,2);

MMQ(1,:) = DATA(1:first_end-first_start+1,2);
MMQ(2:4,:) = DATA(first_start:first_end,6:8)'; %m/s^2
MMQ(5:7,:) = DATA(first_start:first_end,3:5)'; %rad/s

MMQ(3,:)=-MMQ(3,:);
MMQ(4,:)=-MMQ(4,:);

MMQ(6,:)=-MMQ(6,:);
MMQ(7,:)=-MMQ(7,:);

% first_start=1;
% first_end=50001;
% first_step=first_end-first_start;
% gyro_offsetx=sum(MMQ(5,first_start:first_end))/first_step;
% gyro_offsety=sum(MMQ(6,first_start:first_end))/first_step;
% gyro_offsetz=sum(MMQ(7,first_start:first_end))/first_step;
% 
% acce_offsetx=sum(MMQ(2,first_start:first_end))/first_step;
% acce_offsety=sum(MMQ(3,first_start:first_end))/first_step;
% acce_offsetz=sum(MMQ(4,first_start:first_end))/first_step;


n=length(time);        	

tmp1=distance(1);
tmp2=time(1);

v=zeros();
vv=zeros();
j=1;

start_n=21;
step_n=20;
delta_t=0.005*step_n;
Ds=zeros();

for i=start_n:step_n:n
    ds1=distance(i-step_n+5)-distance(i-step_n);
    dv1=ds1/(5*0.005);
    ds2=distance(i-step_n+10)-distance(i-step_n+5);
    dv2=ds2/(5*0.005);
    ds3=distance(i-step_n+15)-distance(i-step_n+10);
    dv3=ds3/(5*0.005);
    ds4=distance(i-step_n+20)-distance(i-step_n+15);
    dv4=ds4/(5*0.005);
    
    ds=distance(i)-tmp1;
    tmp1=distance(i);
    dt=time(i)-tmp2;
    tmp2=time(i);
    Ds(j+1)=Ds(j)+ds;
%     Ds(j+1)=Ds(j)+ds;
    v(j)=(dv1+dv2+dv3+dv4)/4;
    vv(j)=ds/dt;
    j=j+1;
end
% deep=10;
% buf_v=zeros(deep,1);
% smooth_v=v;
% for i=1:j-1
%     for k=deep:-1:2
%         buf_v(k)=buf_v(k-1);
%     end
%     buf_v(1)=v(i);
%     smooth_v(i)=sum(buf_v)/deep;
%     
% end
GPS(1,:) = [1:j-1]*delta_t;
GPS(2,:) =  40.006403*pi/180;
GPS(3,:) = 116.336532*pi/180;
GPS(4,:) = 50.0;

% GPS(5,:) = smooth_v;
GPS(5,:) = v;
% GPS(6,:) = diff(Ds(1:end))/delta_t;
% GPS(7,:) = smooth_v;
GPS(6,:) = Ds(1:end-1);
% GPS(7,:) = diff(Ds(1:end-1))*10;
save('coreTurn.mat','MMQ','GPS')		


