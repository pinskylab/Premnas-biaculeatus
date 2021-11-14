

cis = data.frame(i=1:500, lower=NA)
cis$lower[i]=De*i/qchisq(0.025, df=i)
  
for(i in 1:500) {
  cis=c(De*i/qchisq(0.025, df=i))  #Need to add df= otherwise get an error message
  print(cis)
}

#df=3 is 815.1128, seems closest to lower CI of 779 from Ne estimator
#df=4 is 484.1474

for(i in 1:500) {
  cis=c(De*i/qchisq(0.025, df=i))
}

cis 

#When I try to append the cis to the data frame, it only prints the value for i=500  
#print(cis) gives me everything, but doesn't list the i next to it
cis_df <- data.frame(i=1:500, lower=cis)
cis_df

