#this script is to simulate different uni- and multivariate datasets and test statistical methods, e.g. method by Lebl et al. 2013
# using cross-correlation maps and Poisson regression and genetic algorithms
library(lattice)

rm(list = ls());

#create univariate time-series of length 10*52, sinusoidal with wavelength of 52 weeks, amplitude of 1
t=seq(0,520,1);
Sinuswave_expl<-ts(1*sin(2*pi*1/52*t))
plot(t/52,Sinuswave_expl)

#create another univariate time-series of length 10*52, sinusoidal with wavelength of 52 weeks, amplitude of 1
t=seq(0,520,1);
X2<-ts(1*sin(2*pi*1/52*t))
plot(t/52,X2)

#create a random noise explanatory variable with mean 0 and sd 1 with same length
t=seq(0,520,1);
Random_expl<-ts(rnorm(length(t)));
plot(t/52,Random_expl)

#create a random noise response variable with mean 0 and sd 1 with same length
t=seq(0,520);
Random_resp<-ts(rnorm(length(t)));
plot(t/52,Random_resp)

#create response variable perfectly related to X1 with a time lag of 4 weeks
Sinuswave_resp_lag=Sinuswave_expl[(length(Sinuswave_expl)-3):length(Sinuswave_expl)];
Sinuswave_resp_lag<-c(Sinuswave_resp_lag,Sinuswave_expl[1:(length(Sinuswave_expl)-4)]);
Sinuswave_resp_lag<-ts(Sinuswave_resp_lag);

#now let's see the cross-correlation function between Y2 and X1
Sinuswave_resp_lag_Sinuswave_expl_ccf<-ccf(Sinuswave_expl,Sinuswave_resp_lag,lag.max=40)
plot(Sinuswave_resp_lag_Sinuswave_expl_ccf)

#now let's try the opposite
Sinuswave_expl_Sinuswave_resp_lag_ccf<-ccf(Sinuswave_resp_lag,Sinuswave_expl,lag.max=40)
plot(Sinuswave_expl_Sinuswave_resp_lag_ccf)

#now let's try what we get between Y2 and random explanatory X3
Sinuswave_resp_lag_Random_expl_ccf<-ccf(Random_expl,Sinuswave_resp_lag,lag.max=40)
plot(Sinuswave_resp_lag_Random_expl_ccf)

#now let's try what we get between random Y1 and explanatory X1
Random_resp_Sinuswave_expl_ccf<-ccf(Sinuswave_expl,Random_resp,lag.max=40)
plot(Random_resp_Sinuswave_expl_ccf)

#now let's try what we get between random Y1 and random X3
Random_resp_Random_expl_ccf<-ccf(Random_expl,Random_resp,lag.max=40)
plot(Random_resp_Random_expl_ccf) #note that there are significant peaks (by random)

#let's make a new explanatory variable with sum of sine wave and random noise w 1/5th weight
Mixed_expl=Sinuswave_expl+Random_expl/5
plot(Mixed_expl)

#let's make a new response variable with sum of sine wave and random noise w 1/5th weight
Mixed_resp_lag=Sinuswave_resp_lag+Random_resp/5
plot(Mixed_resp_lag)

#now let's try what we get between sine wave Y2 and mixed X4 # still pretty good
Sinuswave_resp_lag_Mixed_expl_ccf<-ccf(Mixed_expl,Sinuswave_resp_lag,lag.max=40)
plot(Sinuswave_resp_lag_Mixed_expl_ccf)

#now let's try what we get between mixed Y3 and sine wave X1 #still pretty good
Mixed_resp_lag_Sinuswave_expl_ccf<-ccf(Sinuswave_expl,Mixed_resp_lag,lag.max=40)
plot(Mixed_resp_lag_Sinuswave_expl_ccf)

#now let's try what we get between mixed Y3 and mixed X4 #still pretty good #still pretty good
Y3_X4_ccf<-ccf(X4,Y3,lag.max=40)
plot(Y3_X4_ccf) #strong signal

#now what if I difference both Y2 and X1
diff_X1<-diff(X1)
plot(diff_X1)
diff_Y2<-diff(Y2)
lines(diff_Y2,col="red")

#let's look at ccf between differenced Y2 and X1
diff_Y2_diff_X1_ccf<-ccf(diff_X1,diff_Y2,lag.max=40)
plot(diff_Y2_diff_X1_ccf)   #still the same, working good

#now what if I difference both Y3 and X4
diff_X4<-diff(X4)
plot(diff_X4)
diff_Y3<-diff(Y3)
lines(diff_Y3,col="red")

#let's look at ccf between differenced Y3 and X4
diff_Y3_diff_X4_ccf<-ccf(diff_X4,diff_Y3,lag.max=40)
plot(diff_Y3_diff_X4_ccf)   #more noisy, hard to pick out exact lag #noise destroys information in cross-correlation function

#let's try to do cross-correlation map between Y2 and X1
max.lag=20;
CCM_matrix<-NULL;
for (i in (max.lag+1):length(X1))
for (j in 0:max.lag)
for (k in j:max.lag)
CCM_matrix<-rbind(CCM_matrix,c(i,j,k,mean(X1[(i-k):(i-j)],na.rm=TRUE)));

CCM_corr_matrix<-NULL;
for (j in 0:max.lag)
for (k in 0:max.lag)
{
	relevant_rows<-which(CCM_matrix[,2]==j&CCM_matrix[,3]==k);
	if (length(relevant_rows)>0)
	{
	estimate_Pearson<-cor.test(CCM_matrix[relevant_rows,4],Y2[CCM_matrix[relevant_rows,1]],method="pearson")$estimate;
	estimate_Spearman<-cor.test(CCM_matrix[relevant_rows,4],Y2[CCM_matrix[relevant_rows,1]],method="spearman")$estimate;	
	}
	if (length(relevant_rows)==0) {estimate_Pearson=NA;estimate_Spearman=NA};
CCM_corr_matrix<-rbind(CCM_corr_matrix,c(j,k,estimate_Pearson,estimate_Spearman));
}

rgb.palette <- colorRampPalette(c("blue", "cyan", "green","yellow", "red" , "darkred"), space = "rgb")

#let's visualize CCM for Pearson correlation

Pearson_corr_matrix<-matrix(CCM_corr_matrix[,3], nrow=21, ncol=21, dimnames=list(c(1:21), c(1:21)))
levelplot(Pearson_corr_matrix, xlab="", ylab="", main="Pearson Y2~X1", col.regions=rgb.palette(100), cuts=40)

#let's visualize CCM for Spearman correlation

Spearman_corr_matrix<-matrix(CCM_corr_matrix[,4], nrow=21, ncol=21, dimnames=list(c(1:21), c(1:21)))
levelplot(Spearman_corr_matrix, xlab="", ylab="", main="Spearman Y2~X1", col.regions=rgb.palette(100), cuts=40)

#let's compare Pearson's to Spearman's correlation #almost the same in this case
plot(CCM_corr_matrix[,3],CCM_corr_matrix[,4],xlab="Pearson's",ylab="Spearman's")
lm_comparison<-lm(CCM_corr_matrix[,3]~CCM_corr_matrix[,4])

#maximum of the corrplot (Spearmans) correctly identifies the lag 4,4 as the highest correlation
CCM_corr_matrix[which(CCM_corr_matrix[,4]==max(CCM_corr_matrix[,4],na.rm=T)),]

#maximum of the corrplot (Pearsons) correctly identifies the lag 4,4 as the highest correlation
CCM_corr_matrix[which(CCM_corr_matrix[,3]==max(CCM_corr_matrix[,3],na.rm=T)),]

#let's check if the diagonal of the Pearson's correlation matches with the ccf
diag_Pearson<-Pearson_corr_matrix[seq(1,length(Pearson_corr_matrix),22)];
ccf_estimates<-Y2_X1_ccf$acf[41:21]
plot(diag_Pearson,ccf_estimates,xlab="Pearson's",ylab="ccf")

#let's try the CCM with the noisy dataset Y3 and X4
max.lag=20;
CCM_matrix<-NULL;
for (i in (max.lag+1):length(X4))
for (j in 0:max.lag)
for (k in j:max.lag)
CCM_matrix<-rbind(CCM_matrix,c(i,j,k,mean(X4[(i-k):(i-j)],na.rm=TRUE)));

CCM_corr_matrix<-NULL;
for (j in 0:max.lag)
for (k in 0:max.lag)
{
	relevant_rows<-which(CCM_matrix[,2]==j&CCM_matrix[,3]==k);
	if (length(relevant_rows)>0)
	{
	estimate_Pearson<-cor.test(CCM_matrix[relevant_rows,4],Y3[CCM_matrix[relevant_rows,1]],method="pearson")$estimate;
	estimate_Spearman<-cor.test(CCM_matrix[relevant_rows,4],Y3[CCM_matrix[relevant_rows,1]],method="spearman")$estimate;	
	}
	if (length(relevant_rows)==0) {estimate_Pearson=NA;estimate_Spearman=NA};
CCM_corr_matrix<-rbind(CCM_corr_matrix,c(j,k,estimate_Pearson,estimate_Spearman));
}

#let's visualize CCM for Pearson correlation

Pearson_corr_matrix<-matrix(CCM_corr_matrix[,3], nrow=21, ncol=21, dimnames=list(c(1:21), c(1:21)))
levelplot(Pearson_corr_matrix, xlab="", ylab="", main="Pearson Y3~X4", col.regions=rgb.palette(100), cuts=40)

#let's visualize CCM for Spearman correlation

Spearman_corr_matrix<-matrix(CCM_corr_matrix[,4], nrow=21, ncol=21, dimnames=list(c(1:21), c(1:21)))
levelplot(Spearman_corr_matrix, xlab="", ylab="", main="Spearman Y3~X4", col.regions=rgb.palette(100), cuts=40)

#let's compare Pearson's to Spearman's correlation #almost the same in this case
plot(CCM_corr_matrix[,3],CCM_corr_matrix[,4],xlab="Pearson's",ylab="Spearman's")
lm_comparison<-lm(CCM_corr_matrix[,3]~CCM_corr_matrix[,4])

#let's check if the diagonal of the Pearson's correlation matches with the ccf
diag_Pearson<-Pearson_corr_matrix[seq(1,length(Pearson_corr_matrix),22)];
ccf_estimates<-Y3_X4_ccf$acf[41:21]
plot(diag_Pearson,ccf_estimates,xlab="Pearson's",ylab="ccf")

#maximum of the corrplot (Spearmans) correctly identifies the lag 4,4 as the highest correlation
CCM_corr_matrix[which(CCM_corr_matrix[,4]==max(CCM_corr_matrix[,4],na.rm=T)),]

#maximum of the corrplot (Pearsons) correctly identifies the lag 4,4 as the highest correlation
CCM_corr_matrix[which(CCM_corr_matrix[,3]==max(CCM_corr_matrix[,3],na.rm=T)),]

#check if the ccf is doing what it's supposed to be doing
cor.test(Y3[1:length(Y3)],X4[1:(length(X4)-1)])  #ccf works with lags of x[t+k] and y[t]