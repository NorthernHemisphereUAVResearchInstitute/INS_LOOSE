function [nRn,nRe,gn,Wie,Wep] = EarthCalc(lat,he,V)

  if nargin<3,error('insufficient number of input arguments'),end
    e  =0.003352813;
    We =0.00007292115;
    Re =6378137.0;
    
    nRn=Re*(1-2*e+3*e*sin(lat)*sin(lat));
    nRe=Re*(1+e*sin(lat)*sin(lat));
    g=normgravity(lat,he);
    gn = [0;0;g];

    Wie = [We*cos(lat);0;-We*sin(lat)];
    Wep = [V(2)/(nRe+he);-V(1)/(nRn+he);-V(2)*tan(lat)/(nRe+he)];
