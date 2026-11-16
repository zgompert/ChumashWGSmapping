## calculation of demogrpahic parameters
## from momments SI model
theta_si<-665146.2194320768
T_si<-0.5750388832301452

## meidan and 95% HDPIs from Victor's tree, based on 40,000 trees from posterior
## divergence time in millions of years between HF and GR8
Time<-c(5300000,2256000,9580000)

## Ne under SI, theta = 4 Ne u
## T_si = 2 Ne generation
## thus Time = T_si * 2Ne
## Ne = Time/(2*T_si)

Ne_si<-Time/(2*T_si)
Ne_si
#[1] 4608384 1961606 8329871

mu<-theta_si/(4*Ne_si)
#[1] 0.03608348 0.08477060 0.01996268

## now for ancient migration
theta_am<-236974.3856507005
nu1<-1.5170759930354294
nu2<-3.598688723355489
T1<-3.989565807508886
T2<-0.024923114593192987
M12<-0.035697270880724004## from 2 into 1
M21<-0.07219446440497963## from 1 into 2

Ne_am<-theta_am/(4*mu)
#[1] 1364464.4  580798.4 2466333.8

## split time in generations
T1_am<-(T1+T2) * 2 * Ne_am
#[1] 13182361  5611209 23827739

## time migration halts
T2_am<-T2 * 2 * Ne_am
#[1]  81839.93  34836.02 147929.53

## Ne in pops 1 and 2 after split
Ne_am_1<-nu1*Ne_am
#[1] 2490808 1060238 4502253
Ne_am_2<-nu2*Ne_am
#[1] 5908500  2515014 10679892

## migration in number of migrants
M12_am<-M12/2 * nu2
#[1] 0.06423168
M21_am<-M21/2 * nu1
#[1] 0.05476224


