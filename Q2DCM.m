function Cnb=Q2DCM(Qnb)

Cnb(1,1)=Qnb(1)*Qnb(1)+Qnb(2)*Qnb(2)-Qnb(3)*Qnb(3)-Qnb(4)*Qnb(4);
Cnb(1,2)=2*(Qnb(2)*Qnb(3)-Qnb(1)*Qnb(4));
Cnb(1,3)=2*(Qnb(2)*Qnb(4)+Qnb(1)*Qnb(3));
Cnb(2,1)=2*(Qnb(2)*Qnb(3)+Qnb(1)*Qnb(4));
Cnb(2,2)=Qnb(1)*Qnb(1)-Qnb(2)*Qnb(2)+Qnb(3)*Qnb(3)-Qnb(4)*Qnb(4);
Cnb(2,3)=2*(Qnb(3)*Qnb(4)-Qnb(1)*Qnb(2));
Cnb(3,1)=2*(Qnb(2)*Qnb(4)-Qnb(1)*Qnb(3));
Cnb(3,2)=2*(Qnb(3)*Qnb(4)+Qnb(1)*Qnb(2));
Cnb(3,3)=Qnb(1)*Qnb(1)-Qnb(2)*Qnb(2)-Qnb(3)*Qnb(3)+Qnb(4)*Qnb(4);
