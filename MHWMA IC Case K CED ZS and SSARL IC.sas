proc iml;
**********************************************************************;
*/Zero-state (tau = 1) and steady-state (tau>1)/*;
tau=10;
*Number of simulations;
sim = 10000;
crlvec = j(sim,1,.);
p=3;
*Size of the Phase II sample;
n = 1;

*Design parameters of ;

 mu0 = {0,0,0};
 Sig0 = {1 0 0,0 1 0, 0 0 1};
 *Smoothing parameter;
 lambda1 = 0.25;
 l1 = 1-lambda1;

*Shift;
*d=j(1,3,0.1443376);
d=j(1,p,0.00);
d1=j(1,p,0.00);
delta0 = d*Sig0;
delta1= d1*Sig0;
meanshift0  = mu0+delta0`;
meanshift1  = mu0+delta1`;
***************Control limit************;
h=12.064;

do k=1 to sim; 
count=0;

signal=0;
count2=0;
vec={};

zi_1 = mu0;
zi_2=mu0;
vec={};

df=5;
v=1;
wp=2;

do i=1 to 1000000000 until (signal=1); 
	        yi0 = j(n,1,.);
			yi1 = j(n,1,.);
	        * Generating a Phase II sample;
	        * Generating observations from the Normal distribution;
            yi0= randnormal(n,meanshift0,Sig0);
			yi1= randnormal(n,meanshift1,Sig0);
			if i<tau then yi = yi0; 
	        else yi= yi1;
            count=count+1;
            x_vecb = mean(yi);
			
            vec=vec//x_vecb`;
           **Obtaining MEHWMA the charting statistics;
			zi = lambda1*x_vecb` + (1-lambda1)*zi_1;          
            if i=1 then do; Sig0e=lambda1**2*Sig0/n;
                           	goto test;
							end;
			Sig0e=lambda1**2*Sig0/n+(l1)**2*Sig0/(n*(i-1));
            
            z2b=mean(vec);
	        zi_1 = z2b`;

            test: ti2=((zi-mu0)`*inv(Sig0e))*(zi-mu0); 			
            *Comparing the plotting statistics to the control limits;

            if ((ti2>=h)&(i>tau)) then signal=1;
			else signal=0;
            
            crlvec[k,1]=count-tau;
			end;
   end;

results=crlvec;
print 'HWMA scheme';
print sim n tau lambda1 h;
create runlength_prec_normal from results[colname={d000}];
append from results;
quit;
data  runlength_prec_normal;
set work.runlength_prec_normal;

proc univariate data= runlength_prec_normal noprint;
*var CARLo CSDRLo; 
var d000;		   
histogram;
inset mean std p5 q1 median q3 p95 / format = 10.2;
run;/*
proc print data = pctls2;
run;
