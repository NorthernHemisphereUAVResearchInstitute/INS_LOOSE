function [LI,I] = LIKalman(LI,I,GPS)

global S;

        F = zeros(15,15);
        F(2,3) = I.Wie(1)+I.Wep(1);
        F(1,3) = -(I.Wie(2)+I.Wep(2));
        F(1,2) = I.Wie(3)+I.Wep(3);
        F(2,1) = -(I.Wie(3)+I.Wep(3));
        F(3,1) = I.Wie(2)+I.Wep(2);
        F(3,2) = -(I.Wie(1)+I.Wep(1));
        
        F(1,5) = 1/(I.nRe+LI.height);
        F(2,4) = -1/(I.nRn+LI.height);
        F(3,5) = -tan(LI.latitude)/(I.nRe+LI.height);
        
        F(1,7) = -S.We*sin(LI.latitude);
        F(1,9) =-LI.vn(2)/((I.nRn+LI.height)*(I.nRn+LI.height));
        F(2,9) =LI.vn(1)/((I.nRn+LI.height)*(I.nRn+LI.height));
        F(3,7) = -S.We*cos(LI.latitude)-LI.vn(2)/((I.nRe+LI.height)*cos(LI.latitude)*cos(LI.latitude));
        F(3,9) = LI.vn(2)*tan(LI.latitude)/((I.nRn+LI.height)*(I.nRn+LI.height));
        
        F(1:3,10:12) = -I.Cnb;
        
        F(5,3) = -I.dv_f_n(1)/I.dt;
        F(4,3) = I.dv_f_n(2)/I.dt;
        F(4,2) = -I.dv_f_n(3)/I.dt;
        F(5,1) = I.dv_f_n(3)/I.dt;
        F(6,1) = -I.dv_f_n(2)/I.dt;
        F(6,2) = I.dv_f_n(1)/I.dt;
        
        F(4,4) = LI.vn(3)/(I.nRn+LI.height);
        F(4,5) = -2*(S.We*sin(LI.latitude)+LI.vn(2)*tan(LI.latitude)/(I.nRe+LI.height));
        F(4,6) = LI.vn(1)/(I.nRn+LI.height);
        F(5,4) = 2*S.We*sin(LI.latitude)+LI.vn(2)*tan(LI.latitude)/(I.nRe+LI.height);
        F(5,5) = LI.vn(1)*tan(LI.latitude)/(I.nRe+LI.height)+LI.vn(3)/(I.nRe+LI.height);
        F(5,6) = 2*S.We*cos(LI.latitude)+LI.vn(2)/(I.nRe+LI.height);
        F(6,4) = -2*LI.vn(1)/(I.nRn+LI.height);
        F(6,5) = -2*(S.We*cos(LI.latitude)+LI.vn(2)/(I.nRe+LI.height));
        
        F(4,7) = -(2*S.We*cos(LI.latitude)+LI.vn(2)/((I.nRe+LI.height)*cos(LI.latitude)*cos(LI.latitude)))*LI.vn(2);
        F(4,9)=(LI.vn(2)*LI.vn(2)*tan(LI.latitude)-LI.vn(1)*LI.vn(3))/((I.nRn+LI.height)*(I.nRn+LI.height));
        F(5,7) = (2*S.We*cos(LI.latitude)+LI.vn(2)/((I.nRe+LI.height)*cos(LI.latitude)*cos(LI.latitude)))*LI.vn(1)-2*S.We*LI.vn(3)*sin(LI.latitude);
        F(5,9)=-LI.vn(2)*(LI.vn(1)*tan(LI.latitude)+LI.vn(3))/((I.nRn+LI.height)*(I.nRn+LI.height));
        F(6,7) = 2*LI.vn(2)*S.We*sin(LI.latitude);
        F(6,9) = (LI.vn(1)*LI.vn(1)+LI.vn(2)*LI.vn(2))/((I.nRn+LI.height)*(I.nRn+LI.height));
        
        F(4:6,13:15) = I.Cnb;
        
        F(7,4) = 1/(I.nRn+LI.height);
        F(7,9) = -LI.vn(1)/((I.nRn+LI.height)*(I.nRn+LI.height));
        F(8,5) = 1/((I.nRe+LI.height)*cos(LI.latitude));
        F(9,6) = -1;
        
        F(8,7) = LI.vn(2)*tan(LI.latitude)/((I.nRe+LI.height)*cos(LI.latitude));
        F(8,9) = -LI.vn(2)/((I.nRe+LI.height)*(I.nRe+LI.height)*cos(LI.latitude));
        
        LI.FI = F*LI.dt + eye(15,15);
        
        %速度观测矩阵、观测量
%         zeta=(I.Wie+I.Wep-LI.dx(1:3))*LI.dt;
%         Cn=[1,zeta(3)*0.5,-zeta(2)*0.5;...
%             -zeta(3)*0.5,1,zeta(1)*0.5;...
%             zeta(2)*0.5,-zeta(1)*0.5,1];
        Vned=I.Cnb*[GPS(5,LI.indexGPS);0;0];
%         Vned(1)=Vned(1)*GPS(5,LI.indexGPS)/sqrt(Vned(1)*Vned(1)+Vned(2)*Vned(2)+Vned(3)*Vned(3));
%         Vned(2)=Vned(2)*GPS(5,LI.indexGPS)/sqrt(Vned(1)*Vned(1)+Vned(2)*Vned(2)+Vned(3)*Vned(3));
%         Vned(3)=Vned(3)*GPS(5,LI.indexGPS)/sqrt(Vned(1)*Vned(1)+Vned(2)*Vned(2)+Vned(3)*Vned(3));
%         Vned(1)=0;
%         Vned(2)=0;
%         Vned(3)=0;
        LI.H = zeros(3,15);
        LI.H(1:3,4:6) = eye(3,3);
        
        Z = [LI.vn(1)-Vned(1);LI.vn(2)-Vned(2);LI.vn(3)-Vned(3)];
        LI.R = LI.R0(1:3,1:3);
        
         %位置观测矩阵、观测量
%         LGH=I.Cnb*[GPS(6,LI.indexGPS);0;0];
%         LI.H = zeros(3,15);
%         LI.H(1:3,7:9) = eye(3,3);
%         
%         Z = [LI.latitude-LGH(1);LI.longitude-LGH(2);LI.height-LGH(3)];
%         LI.R = LI.R0(1:3,1:3);
        
        %速度、位置观测矩阵、观测量
%         Vned=I.Cnb*GPS(5:7,LI.indexGPS);
%         LGH=I.Cnb*[GPS(8,LI.indexGPS);0;0];
%         LI.H = zeros(6,15);
%         LI.H(1:6,4:9) = eye(6,6);
%         
%         Z = [LI.vn(1)-Vned(1);LI.vn(2)-Vned(2);LI.vn(3)-Vned(3);...
%              LI.latitude-LGH(1);LI.longitude-LGH(2);LI.height-LGH(3)];
%         LI.R = LI.R0(1:6,1:6);

        Xkk_1 = LI.FI*LI.dx;
        Xkk_1 = Xkk_1 - LI.Tdx;
        Pkk_1 = LI.FI*LI.P*LI.FI'+LI.Q;
        LI.K = Pkk_1*LI.H'*inv(LI.H*Pkk_1*LI.H'+LI.R);
        LI.P = (eye(15,15)-LI.K*LI.H)*Pkk_1;
        LI.dx = Xkk_1 + LI.K*(Z-LI.H*Xkk_1);
        %update
        I.Qnb = qmul(rv2q(LI.dx(1:3,:)),I.Qnb);
        I.Qnb=I.Qnb/norm(I.Qnb);%四元数归一化
        I.Cnb=Q2DCM(I.Qnb);	
        LI.vn = LI.vn - LI.dx(4:6);
        LI.latitude = LI.latitude - LI.dx(7);
        LI.longitude = LI.longitude - LI.dx(8);
        LI.height = LI.height - LI.dx(9);
        LI.Gyro_Offset = LI.Gyro_Offset + LI.dx(10:12);
        LI.Acce_Offset = LI.Acce_Offset + LI.dx(13:15);
        I.Cen=[-sin(LI.latitude)*cos(LI.longitude), -sin(LI.latitude)*sin(LI.longitude),  cos(LI.latitude);
            -sin(LI.longitude),                cos(LI.longitude),                 0;
            -cos(LI.latitude)*cos(LI.longitude),  -cos(LI.latitude)*sin(LI.longitude),  -sin(LI.latitude)];
        
        LI.alllen=GPS(6,LI.indexGPS);
        %save data
        LI.save_dx(:,LI.indexGPS) = LI.dx;
        LI.save_Gyro_Offset(:,LI.indexGPS) =LI.Gyro_Offset;
        LI.save_Acce_Offset(:,LI.indexGPS) =LI.Acce_Offset;
        if(3 == length(Z))
            Z = [Z;0;0;0;0;0];
        else
            Z = [Z;1];
        end
        LI.save_Z(:,LI.indexGPS) = Z;
        LI.save_P(:,LI.indexGPS) = diag(LI.P);
        %set zero

        LI.dx = zeros(15,1);
        LI.Tdx = zeros(15,1);


        