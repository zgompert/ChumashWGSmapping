## calculation of demogrpahic parameters
## from momments SI model
theta_si<-813566.3121126157
T_si<-0.5440668702238484

## meidan and 95% HDPIs from Victor's tree, based on 40,000 trees from posterior
## divergence time in millions of years between HF and GR8
Time<-c(5300000,2256000,9580000)

## Ne under SI, theta = 4 Ne u
## T_si = 2 Ne generation
## thus Time = T_si * 2Ne
## Ne = Time/(2*T_si)

Ne_si<-Time/(2*T_si)
Ne_si
#[1] 4870725 2073275 8804065

mu<-theta_si/(4*Ne_si)
#[1] 0.04175797 0.09810161 0.02310201

## now for ancient migration
theta_am<-196938.52144932494
nu1<-2.15393538297051
nu2<-4.972966519335093 
T1<-6.3044225117597055
T2<-0.02464538523985764
M12<-0.027927807712651236## from 2 into 1
M21<-0.05173863141434871## from 1 into 2

Ne_am<-theta_am/(4*mu)
#[1] 1179047.5  501873.8 2131184.0

## split time in generations
T1_am<-(T1+T2) * 2 * Ne_am
#[1] 14924544  6352787 26976817

## time migration halts
T2_am<-T2 * 2 * Ne_am
#[1] 58116.16  24737.75 105047.70

## Ne in pops 1 and 2 after split
Ne_am_1<-nu1*Ne_am
#[1] 2539592 1081004 4590433
Ne_am_2<-nu2*Ne_am
#[1] 5863364  2495802 10598307

## migration in number of migrants
M12_am<-M12/2 * nu2
#[1] 0.06944203
M21_am<-M21/2 * nu1
#[1] 0.05572083


