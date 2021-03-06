function [LI,I] = LIInit(inittime)

global S;
LI.longitude = 116.336532*S.PI/180;
LI.latitude = 40.006403*S.PI/180;
LI.height = 50.0;
LI.alllen = 0;

LI.vn(1,1) = 0;
LI.vn(2,1) = 0;
LI.vn(3,1) = 0;
LI.sita = 4.670000000000000*S.Con_d2r;
LI.gama = 0.010000000000000*S.Con_d2r;
LI.psi = (2.649200000000000e+02 + 5)*S.Con_d2r;   % 磁偏角修正：大于180度要减5度，小于180度要加5度
LI.psi_GPS_flag = 0;
% LI.Gyro_Offset = [1.364698593971882e-05/3600;-5.239851817036305e-05/3600;4.885449968999365e-05/3600];%[30*S.PI/180/3600;30*S.PI/180/3600;30*S.PI/180/3600];%[0.00097201999999999998;-0.0024828200000000002;0.0014251200000000000];
LI.Gyro_Offset = [0;0;0];
LI.Acce_Offset = [0;0;0];%[-0.10568246599;-9.736858616700014+9.7803267715;0.546513824656927][1e-3*9.8;1e-3*9.8;1e-3*9.8];%[0.031241300000000000;-0.089116100000000004;-0.0039435599999999996];

LI.dt = 0.1;%kalman滤波器滤波周期
LI.time = inittime;

I.oaSet = [0;0;0];%人为加入零偏
I.ogSet = [0;0;0];%[2*1.4544e-004; 2*1.4544e-004; 2*1.4544e-004];
I.Wie = [S.We*cos(LI.latitude);0;-S.We*sin(LI.latitude)];
I.Cen=[-sin(LI.latitude)*cos(LI.longitude), -sin(LI.latitude)*sin(LI.longitude),  cos(LI.latitude);
    -sin(LI.longitude),                cos(LI.longitude),                 0;
    -cos(LI.latitude)*cos(LI.longitude),  -cos(LI.latitude)*sin(LI.longitude),  -sin(LI.latitude)];

I.Cnb(1,1)=cos(LI.sita)*cos(LI.psi);
I.Cnb(1,2)=-cos(LI.gama)*sin(LI.psi)+sin(LI.gama)*sin(LI.sita)*cos(LI.psi);
I.Cnb(1,3)=sin(LI.gama)*sin(LI.psi)+cos(LI.gama)*sin(LI.sita)*cos(LI.psi);
I.Cnb(2,1)=cos(LI.sita)*sin(LI.psi);
I.Cnb(2,2)=cos(LI.gama)*cos(LI.psi)+sin(LI.gama)*sin(LI.sita)*sin(LI.psi);
I.Cnb(2,3)=-sin(LI.gama)*cos(LI.psi)+cos(LI.gama)*sin(LI.sita)*sin(LI.psi);
I.Cnb(3,1)=-sin(LI.sita);
I.Cnb(3,2)=sin(LI.gama)*cos(LI.sita);
I.Cnb(3,3)=cos(LI.gama)*cos(LI.sita);
I.Qnb=Cnb2qnb(I.Cnb);
I.Qnb=I.Qnb/norm(I.Qnb);%四元数归一化

I.Dvn=zeros(3,1);

I.fs = 200;%采样频率
I.ts=1/I.fs;%采样时间
I.dt = 2*I.ts;
I.dv_f_n = zeros(3,1);
%根据纬度和高度对地球参数进行修正
[I.nRn,I.nRe,I.gn,I.Wie,I.Wep] = EarthCalc(LI.latitude,LI.height,LI.vn);
%初始化kalman滤波参数
LI.FI = eye(15,15);
%LI.G = eye(15,15);
%% 原始参数
% LI.R = diag([(0.001*S.PI/180/3600)*(0.001*S.PI/180/3600), (0.001*S.PI/180/3600)*(0.001*S.PI/180/3600), ...
%     0.01*0.01, 0.1*0.1, 0.1*0.1, 0.5*0.5]);%观测噪声方差阵——观测量为位置
% LI.P = diag([0.2*S.Con_d2r,0.2*S.Con_d2r,0.2*S.Con_d2r,0.1,0.1,0.1,0.01*S.PI/180/3600,...
%     0.01*S.PI/180/3600,0.1,1*S.Con_dh2rs,1*S.Con_dh2rs,1*S.Con_dh2rs,0.01,0.01,0.01]);
% LI.Q = diag([(60.0/3600*S.Con_d2r)*LI.dt, (60.0/3600*S.Con_d2r)*LI.dt, ...
%     (60.0/3600*S.Con_d2r)*LI.dt, 0.002*LI.dt, 0.002*LI.dt,0.002*LI.dt, ...
%     0, 0, 0, 0, 0, 0, 0, 0, 0]);
% %% 包超的参数
% LI.R = diag([(0.01*S.PI/180/3600)*(0.01*S.PI/180/3600), (0.01*S.PI/180/3600)*(0.01*S.PI/180/3600), ...
%     0.2*0.2, 0.1*0.1, 0.1*0.1, 0.5*0.5]);%观测噪声方差阵——观测量为位置（输入）
% 
% LI.P = diag([0.1^2,0.1^2,0.1^2,2^2,2^2,2^2,(10*S.PI/180/3600)^2,...
%     (10*S.PI/180/3600)^2,3^2,(0.5*S.Con_dh2rs)^2,(0.5*S.Con_dh2rs)^2,(0.5*S.Con_dh2rs)^2,0.005^2,0.005^2,0.005^2]);%状态协方差矩阵初始值（衡量初始状态下噪声大小）
% 
% LI.Q = diag([(120.0/3600*S.Con_d2r)^2*LI.dt, (120.0/3600*S.Con_d2r)^2*LI.dt, ...
%     (750/3600*S.Con_d2r)^2*LI.dt, 0.004^2*LI.dt, 0.004^2*LI.dt,0.004^2*LI.dt, ...
%     0, 0, 0, 0, 0, 0, 0, 0, 0]); %状态噪声方差阵（输入）
%% 低精度IMU参数
% LI.R0 = diag([(0.1*S.PI/180/3600)*(0.1*S.PI/180/3600), (0.1*S.PI/180/3600)*(0.1*S.PI/180/3600), ...
%     5*5, 0.1*0.1, 0.1*0.1, 0.1*0.1, (0.5*S.PI/180.0)*(0.5*S.PI/180.0)]);%观测噪声方差阵——观测量为位置（输入）
LI.R0 = diag([0.5*0.5, 0.5*0.5, 0.5*0.5, ...
    1*1, 1*1, 1*1, (0.5*S.PI/180.0)*(0.5*S.PI/180.0)]);
% LI.R0 = diag([0.1*0.1, 0.1*0.1, ...
%     0.2*0.2, 0.1*0.1, 0.1*0.1, 0.2*0.2, (0.5*S.PI/180.0)*(0.5*S.PI/180.0)]);%观测噪声方差阵——观测量为位置（输入）

LI.P = diag([5^2,5^2,5^2,2^2,2^2,2^2,(10*S.PI/180/3600)^2,...
    (10*S.PI/180/3600)^2,30^2,(0.5*S.Con_dh2rs)^2,(0.5*S.Con_dh2rs)^2,(0.5*S.Con_dh2rs)^2,0.005^2,0.005^2,0.005^2]);%状态协方差矩阵初始值（衡量初始状态下噪声大小）

% LI.Q = diag([(12.0/60*S.Con_d2r)^2, (12.0/60*S.Con_d2r)^2, (12/60*S.Con_d2r)^2,  ...
%     (18/60)^2, (18/60)^2,(18/60)^2, ...
%     0, 0, 0, ...
%     (93/3600*S.Con_d2r)^2, (93/3600*S.Con_d2r)^2, (93/3600*S.Con_d2r)^2,...
%     (40000 * 1.0e-6 * 9.8)^2, (40000 * 1.0e-6 * 9.8)^2, (40000 * 1.0e-6 * 9.8)^2]); %状态噪声方差阵（输入）

% LI.Q = diag([(1.2*S.Con_d2r)^2, (1.2*S.Con_d2r)^2, (1.2*S.Con_d2r)^2,  ...
%     (0.2)^2, (0.2)^2,(0.2)^2, ...
%     0, 0, 0, ...
%     (93/3600*S.Con_d2r)^2, (93/3600*S.Con_d2r)^2, (93/3600*S.Con_d2r)^2,...
%     (10000 * 1.0e-6 * 9.8)^2, (10000 * 1.0e-6 * 9.8)^2, (10000 * 1.0e-6 * 9.8)^2]); %状态噪声方差阵（输入）

LI.Q = diag([(1.2*S.Con_d2r)^2, (1.2*S.Con_d2r)^2, (1.2*S.Con_d2r)^2,  ...
    (0.2)^2, (0.2)^2,(0.2)^2, ...
    0, 0, 0, ...
    0, 0, 0, ...
    0, 0, 0]); %状态噪声方差阵（输入）

% LI.Q = diag([(0.2)^2, (0.2)^2,(0.2)^2,  ...
%     (0.2)^2, (0.2)^2,(0.2)^2, ...
%     0, 0, 0, ...
%     0, 0, 0, ...
%     0, 0, 0]); %状态噪声方差阵（输入）

%% 用于debug的参数
% LI.R = diag([(0.01*S.PI/180/3600)*(0.01*S.PI/180/3600), (0.01*S.PI/180/3600)*(0.01*S.PI/180/3600), ...
%     0.2*0.2, 1*1, 1*1, 5*5]);%观测噪声方差阵——观测量为位置（输入）
% 
% LI.P = diag([0.2*S.Con_d2r,0.2*S.Con_d2r,0.2*S.Con_d2r,0.4,0.4,0.4,0.1*S.PI/180/3600,...
%     0.1*S.PI/180/3600,0.2,2.0*S.Con_dh2rs,2.0*S.Con_dh2rs,2.0*S.Con_dh2rs,0.02,0.02,0.02]);%状态协方差矩阵初始值（衡量初始状态下噪声大小）
% 
% LI.Q = diag([(120.0/3600*S.Con_d2r)*LI.dt, (120.0/3600*S.Con_d2r)*LI.dt, ...
%     (300.0/3600*S.Con_d2r)*LI.dt, 0.004*LI.dt, 0.004*LI.dt,0.004*LI.dt, ...
%     0, 0, 0, 0, 0, 0, 0, 0, 0]); %状态噪声方差阵（输入）

%LI.R = LI.R*10; LI.P = LI.P*10; LI.Q = LI.Q*10;
LI.dx = zeros(15,1);
LI.Tdx = zeros(15,1);
LI.indexGPS = 1;
LI.indexINS = 1;
LI.flagRun = true;