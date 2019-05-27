function q = rv2q(rv)
    F1	=(   2 * 1);		%/* define: Fk=2^k*k! */
    F2	=(F1*2 * 2);
    F3	=(F2*2 * 3);
    F4	=(F3*2 * 4);
    F5	=(F3*2 * 5);
	n2 = rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3);
	if(n2<0.017^2)	%/* 0.017^2 */
		n4=n2*n2;
		c = 1.0 - n2*(1.0/F2) + n4*(1.0/F4);
		f = 0.5 - n2*(1.0/F3) + n4*(1.0/F5);
	else
		n_2 = sqrt(n2);
		c = cos(n_2);
		f = sin(n_2)/n_2*0.5;
    end
	q = [c; f*rv];

%     norm = sqrt(rv'*rv);
%     if norm>1.e-20
%         f = sin(norm/2)/(norm/2);
%     else
%         f = 1;
%     end
%     q = [cos(norm/2); f/2*rv];
