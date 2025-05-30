#!/usr/bin/perl

$out = "post_gemma_GR8_color.txt";
#$out = "post_gemma_color.txt";

system "touch $out\n";

open(R, "> source.R") or die "faile to write R file\n";

$nph = 2;

print R "phm<-matrix(NA,nrow=$nph,ncol=9)\n";

$p = 1;
foreach $ph (@ARGV){
		print R "phl<-vector(\"list\",60)\n";
		foreach $i (1..60){
			$j = $i-1;
			$ph =~ s/ch\d+/ch$j/;
			$in = $ph;
			print R "phl[[$i]]<-read.table(\"$in\",header=TRUE)\n";
		}
		print R "ph<-c(phl[[1]][,2]";
		foreach $i (2..60){
			print R ",phl[[$i]][,2]";
		}
		print R ")\n
                q<-quantile(ph,probs=c(0.5,0.05,0.95))\n
		phm[$p,1:3]<-q\n";
		print R "ph<-c(phl[[1]][,4]";
		foreach $i (2..60){
			print R ",phl[[$i]][,4]";
		}
		print R ")\n
                q<-quantile(ph,probs=c(0.5,0.05,0.95))\n
		phm[$p,4:6]<-q\n";

		print R "ph<-c(phl[[1]][,6]";
		foreach $i (2..60){
			print R ",phl[[$i]][,6]";
		}
		print R ")\n
                q<-quantile(ph,probs=c(0.5,0.05,0.95))\n
		phm[$p,7:9]<-q\n";
		$p++;	
}	
print R "write.table(phm,\"$out\",quote=FALSE,row.names=FALSE,col.names=FALSE)\n";
print R "rm(list=ls())\n";
close(R);
system "R CMD BATCH source.R\n";

