proc iml;
**********************************************************************;
*/Zero-state (tau = 1) and steady-state (tau>1)/*;
tau=100;
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
 lambda = 0.05;
 l1 = 1-lambda;

*Shift;
*d1=j(1,3,0.1443376);
d0=j(1,p,0.00);
d1=j(1,p,0.00);
delta0 = d0*Sig0;
delta1= d1*Sig0;
meanshift0  = mu0+delta0`;
meanshift1  = mu0+delta1`;
***************Control limit************;
h=9.101;

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
			*Obtaining the charting statistics;

	        *comb = xi // yi;
            x_vecb = mean(yi);

           *Plotting statistic;
			zi = lambda*x_vecb` + (1-lambda)*zi_1;
            Cov0e=lambda*(1-(l1)**(2*i))*Sig0/(2-lambda);

            ti2=((zi-mu0)`*inv(Cov0e))*(zi-mu0);
            zi_1 = zi;

	        *Comparing the plotting statistics to the control limits;
  			vec=vec//zi;			

            if ((ti2>=h)&(i>tau)) then signal=1;
			else signal=0;
            
            crlvec[k,1]=count-tau;
			end;
   end;

results=crlvec;
print sim n tau lambda h;
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
