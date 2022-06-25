proc iml;
**********************************************************************;
finalresults = j(13,50,.);
*/Zero-state (tau = 1) and steady-state (tau>1)/*;
do tau=1 to 50;
*Number of simulations;
sim = 10000;
rlvec = j(sim,13,.);
p=3;
*Size of the Phase II sample;
n = 1;

*Design parameters of ;

 mu0 = {0,0,0};
 Sig0 = {1 0 0,0 1 0, 0 0 1};
 *Smoothing parameter;
 lambda1 = 0.1;
 l1 = 1-lambda1;
 lambda2 = 0.05;
 l2 = 1-lambda2;

***************Control limit************;
h=10.635;

*Shift;
d00=j(1,3,0.00);
d0=j(1,3,0.0000000)//j(1,3,0.1443376)//j(1,3,0.2886752)//j(1,3,0.4330129)//j(1,3,0.5773540)//j(1,3,0.7269040)//j(1,3,0.8660400)//j(1,3,1.0103630)//j(1,3,1.1547030)//j(1,3,1.2990384)//j(1,3,1.4433770)//j(1,3,1.5877134)//j(1,3,1.7320510);
do col=1 to 13;
d=d0[col,];
delta = d*Sig0;
delta0=d00*Sig0;
meanshift  = mu0+delta`;
meanshift0  = mu0+delta0`;

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
*/Zero-state (tau = 1) and steady-state (tau>1)/*;

do i=1 to 1000000000 until (signal=1); 
	        yi0 = j(n,1,.);
			yi1 = j(n,1,.);
	        * Generating a Phase II sample;
	        * Generating observations from the Normal distribution;
			yi0= randnormal(n,meanshift0,Sig0);
	        yi1= randnormal(n,meanshift,Sig0);
			if i<tau then yi = yi0; 
	        else yi= yi1;

            count=count+1;
            x_vecb = mean(yi);
			
            vec=vec//x_vecb`;
           **Obtaining MEHWMA the charting statistics;
			zi = lambda1*x_vecb`-lambda2*zi_1 + (1-lambda1+lambda2)*zi_2;          
            if i=1 then do; Sig0e=lambda1**2*Sig0/n;
                           	goto test;
							end;
			Sig0e=(lambda1**2+((l1-(i-2)*lambda2)/(i-1))**2+((l1+lambda2)/(i-1))**2*(i-2))*Sig0/n;
            
            zi_1=x_vecb`;

            z2b=mean(vec);
	        zi_2 = z2b`;

            test: ti2=((zi-mu0)`*inv(Sig0e))*(zi-mu0); 			
            *Comparing the plotting statistics to the control limits;

            if ((ti2>=h)&(i>tau)) then signal=1;
			else signal=0;
            
            rlvec[k,col]=count-tau;
			end;
   end;

end;
results=rlvec;
name1={"0.00", "0.25", "0.50", "0.75", "1.00", "1.25", "1.50", "1.75","2.00", "2.25", "2.50", "2.75", "3.00"};	
results11=mean(rlvec);
results12=std(rlvec);
results21=results11`;
results22=results12`;
*results2=results21||results22;
finalresults[,tau] =results21;
*print sim n tau lambda1 h;
*print results2[rowname=name1][format=10.2];
end;
print 'MEHWMA CED PROPERTIES';
print sim n lambda1 lambda2 h;
name2={"tau=1" "tau=2" "tau=3" "tau=4" "tau=5" "tau=6" "tau=7" "tau=8" "tau=9" "tau=10" "tau=11" "tau=12" "tau=13" "tau=14" "tau=15" "tau=16" "tau=17" "tau=18" "tau=19" "tau=20" "tau=21" "tau=22" "tau=23" "tau=24" "tau=25" "tau=26" "tau=27" "tau=28" "tau=29" "tau=30" "tau=31" "tau=32" "tau=33" "tau=34" "tau=35" "tau=36" "tau=37" "tau=38" "tau=39" "tau=40" "tau=41" "tau=42" "tau=43" "tau=44" "tau=45" "tau=46" "tau=47" "tau=48" "tau=49" "tau=50"};
print finalresults[colname=name2][format=10.2];
quit;
