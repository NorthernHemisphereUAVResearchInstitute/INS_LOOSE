function Euler=DCM2Euler(Cnb,phiold)
PI =3.14159265358979323846 ;
sita=atan(-Cnb(3,1)/sqrt(Cnb(3,2)*Cnb(3,2)+Cnb(3,3)*Cnb(3,3)));
        
if Cnb(3,1)<-0.999
		phi = phiold;
		psi = atan2((Cnb(2,3)-Cnb(1,2)),(Cnb(1,3)+Cnb(2,2)));
elseif Cnb(3,1)>=0.999
		phi = phiold;
		psi = PI+atan2((Cnb(2,3)+Cnb(1,2)),(Cnb(1,3)-Cnb(2,2)));
else
		phi = atan2(Cnb(3,2),Cnb(3,3));
		psi = atan2(Cnb(2,1),Cnb(1,1));
end;
Euler=[sita;phi;psi];
