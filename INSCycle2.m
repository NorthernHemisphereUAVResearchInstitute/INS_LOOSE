function [LI,I] = INSCycle2(LI,I,MMQ)

global S;

    f1 = MMQ(2:4,2*LI.indexINS) + I.oaSet - (LI.Acce_Offset-LI.Tdx(13:15)*I.dt);
    f2 = MMQ(2:4,2*LI.indexINS+1) + I.oaSet - (LI.Acce_Offset-LI.Tdx(13:15)*I.dt);
    w1 = MMQ(5:7,2*LI.indexINS) + I.ogSet - (LI.Gyro_Offset-LI.Tdx(10:12)*I.dt);
    w2 = MMQ(5:7,2*LI.indexINS+1) + I.ogSet - (LI.Gyro_Offset-LI.Tdx(10:12)*I.dt);
    
    drg1=w1*I.ts;
    drg2=w2*I.ts;
    
    dva1=f1*I.ts;
    dva2=f2*I.ts;
    
    drg=drg1+drg2;
    dva=dva1+dva2;
    
    rvet=drg+2/3*cross(drg1,drg2);
    dvb=dva+0.5*cross(drg,dva)+2/3*(cross(dva1,drg2)+cross(drg1,dva2));
    Vmid=LI.vn+I.Dvn;
    Vpre = LI.vn;
    I.Wep = [Vmid(2)/(I.nRe+LI.height);...
            -Vmid(1)/(I.nRn+LI.height);...
            -Vmid(2)*tan(LI.latitude)/(I.nRe+LI.height)];
    zeta=(I.Wie+I.Wep-LI.Tdx(1:3))*I.dt;
    Cn=[1,zeta(3)*0.5,-zeta(2)*0.5;...
        -zeta(3)*0.5,1,zeta(1)*0.5;...
        zeta(2)*0.5,-zeta(1)*0.5,1];
    %速度更新
    I.dv_f_n=Cn*I.Cnb*dvb;
    dv_g_cor=(I.gn-cross(2*I.Wie+I.Wep,Vmid))*I.dt;
    I.Dvn=I.dv_f_n+dv_g_cor;
    LI.vn=LI.vn+I.Dvn-LI.Tdx(4:6)*I.dt;
    %位置更新
    Vmid=0.5*(LI.vn+Vpre);
    I.Wep = [Vmid(2)/(I.nRe+LI.height);...
            -Vmid(1)/(I.nRn+LI.height);...
            -Vmid(2)*tan(LI.latitude)/(I.nRe+LI.height)];
    Cen_cor=[1,I.Wep(3)*I.dt,-I.Wep(2)*I.dt;...
             -I.Wep(3)*I.dt,1,I.Wep(1)*I.dt;...
             I.Wep(2)*I.dt,-I.Wep(1)*I.dt,1];
    I.Cen=Cen_cor*I.Cen;
%     LI.latitude=LI.latitude+Vmid(1)/(I.nRn+LI.height)*I.dt - LI.Tdx(7)*I.dt;
%     LI.longitude=LI.longitude+Vmid(2)/((I.nRe+LI.height)*cos(LI.latitude))*I.dt - LI.Tdx(8)*I.dt;
%     LI.height=LI.height-Vmid(3)*I.dt - LI.Tdx(9)*I.dt;
    
    LI.latitude=LI.latitude+Vmid(1)*I.dt - LI.Tdx(7)*I.dt;
    LI.longitude=LI.longitude+Vmid(2)*I.dt - LI.Tdx(8)*I.dt;
    LI.height=LI.height-Vmid(3)*I.dt - LI.Tdx(9)*I.dt;
    %姿态更新
    Qz=rv2q(-zeta);
    Qr=rv2q(rvet);
    Qn=qmul(I.Qnb,Qr);
    I.Qnb=qmul(Qz,Qn);
    I.Qnb=I.Qnb/norm(I.Qnb);%四元数归一化
    I.Cnb=Q2DCM(I.Qnb);
