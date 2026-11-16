## calculation of demogrpahic parameters
## from momments SI model
theta_si<-872357.1981892756
T_si<-0.6558012695536133

## meidan and 95% HDPIs from Victor's tree, based on 40,000 trees from posterior
## divergence time in millions of years between HF and GR8
Time<-c(5300000,2256000,9580000)

## Ne under SI, theta = 4 Ne u
## T_si = 2 Ne generation
## thus Time = T_si * 2Ne
## Ne = Time/(2*T_si)

Ne_si<-Time/(2*T_si)
Ne_si
#[1] 4040858 1720033 7304042

mu<-theta_si/(4*Ne_si)
#[1] 0.05397103 0.12679365 0.02985871

## now for secondary contact
theta_sc<-483963.77003644564
nu1<-1.20359713036438
nu2<-2.0701761657752753 
T1<-1.357988754540039
T2<-0.7873437492267552
M12<-0.0657494915926104## from 2 into 1
M21<-0.04758511835463534## from 1 into 2

Ne_sc<-theta_sc/(4*mu)
#[1] 2241776  954235 4052115

## split time in generations
T1_sc<-(T1+T2) * 2 * Ne_sc
#[1] 9618708  4094303 17386268

## time migration halts
T2_sc<-T2 * 2 * Ne_sc
#[1] 3530096 1502622 6380815

## Ne in pops 1 and 2 after split
Ne_sc_1<-nu1*Ne_sc
#[1] 2698195 1148515 4877114
Ne_sc_2<-nu2*Ne_sc
#[1] 4640870 1975435 8388592

## migration in number of migrants
M12_sc<-M12/2 * nu2
#[1] 0.06805652
M21_sc<-M21/2 * nu1
#[1] 0.02863666


