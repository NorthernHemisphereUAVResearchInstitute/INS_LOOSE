function g = normgravity(lat,he)

  if nargin<2,error('insufficient number of input arguments'),end
  
    a1 = 9.7803267715;
	a2 = 0.0052790414;
	a3 = 0.0000232718;
	a4 = -0.000003087691089;
	a5 = 0.000000004397731;
	a6 = 0.000000000000721;
	s2 = sin(lat)*sin(lat);
	s4 = s2 * s2;
    
    g= a1 * (1 + a2*s2 + a3*s4) + (a4 + a5*s2)*he + a6 * he * he;
  