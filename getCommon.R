library(data.table)
g<-as.matrix(fread("mlpntest_EXP_tchum.txt",header=FALSE))
p<-apply(g,1,mean)/2
mean(p > 0.05 & p < 0.95)
#[1] 0.0951862
mean(p > 0.025 & p < 0.975)
#[1] 0.1306527
mean(p > 0.01 & p < 0.99)
#[1] 0.1860095

keep<-as.numeric(p > 0.01 & p < 0.99)

write.table(file="KeepSnps.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
