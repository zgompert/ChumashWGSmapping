## calculation of demogrpahic parameters
## from momments SI model
theta_si<-765647.116641579 
T_si<-0.7024078579561003

## meidan and 95% HDPIs from Victor's tree, based on 40,000 trees from posterior
## divergence time in millions of years between HF and GR8
Time<-c(5300000,2256000,9580000)

## Ne under SI, theta = 4 Ne u
## T_si = 2 Ne generation
## thus Time = T_si * 2Ne
## Ne = Time/(2*T_si)

Ne_si<-Time/(2*T_si)
Ne_si
#[1] 3772737 1605905 6819400

mu<-theta_si/(4*Ne_si)
#[1] 0.05073552 0.11919250 0.02806871

## now for ancient migration
theta_am<-254650.2850629205
nu1<-1.9873475456533802
nu2<-3.716254515026489
T1<-4.495927607605403
T2<-0.041616388458101145
M12<-0.03648731748185945## from 2 into 1
M21<-0.03425031825048589## from 1 into 2

Ne_am<-theta_am/(4*mu)
#[1] 1254792.8  534115.6 2268097.2

## split time in generations
T1_am<-(T1+T2) * 2 * Ne_am
#[1] 11387355  4847146 20583182

## time migration halts
T2_am<-T2 * 2 * Ne_am
#[1] 104439.89  44455.92 188780.03 

## Ne in pops 1 and 2 after split
Ne_am_1<-nu1*Ne_am
#[1] 2493709 1061473 4507497
Ne_am_2<-nu2*Ne_am
#[1] 4663130 1984909 8428827

## migration in number of migrants
M12_am<-M12/2 * nu2
#[1]  0.06779808
M21_am<-M21/2 * nu1
#[1] 0.03403364


