r2smw0 = 150;	r2smwlb = 40;	r2smwub = 300;
r2siw0 = 20;	r2siwlb = 6;	r2siwub = 40;
r2sew0 = 30;	r2sewlb = 6;	r2sewub = 40;
mwf    = 0.15;
t1s0   = 225e-3;  	t1slb   = 50e-3;  	t1sub   = 650e-3;
t1l0   = t10;     	t1llb   = 500e-3; 	t1lub   = t10+1;
kls0    = 0;       	klslb    = 0;      	klsub    = 6;       % exchange rate from long T1 to short T1
v_ic = 0.8;
iwf = (1-mwf)*v_ic;
ewf = 1-mwf-iwf;