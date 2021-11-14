##Dispersal Kernels

#From Katrina's Dispersal Variation Repo:

integrate_kernel_sum1 <- function(d, k, theta){ 
  
  disp <-  d*exp(k)*theta*exp(-(exp(k)*d)^theta)/gamma(1/theta)
  return(disp)
  
}

#Gaussian kernel: theta=2, Lapacian kernel: theta=1
#k sets dispersion of the distribution, will use mean sigma=9.01 km
#sigma is the standard deviation of the dispersal kernel

#From Malin's Current Biology paper supplemental methods:

#Probability of a larva settling distance x away from its parent p(x):

#Gaussian kernel: 
#p(x)= (1/sigma*sqrt(2pi))*exp(-(x^2)/2*sigma^2)
#mean dispersal distance is sigma*sqrt(2/pi)

#Lapacian kernel:
#p(x)=(1/sigma*sqrt2)*exp(-sqrt2*x/sigma)
#mean dispersal distance is sigma/sqrt2


#Trying out Gaussian kernel

sigma=5

#Mean Gaussian dispersal distance:

sigma*sqrt(2/3.14)  #7.1987

gaussian_kernel <- function(sigma, x){
  
  px <- (1/sigma*sqrt(2*3.14))*exp(-(x^2)/2*sigma^2)
  return(px)
}

#Want plot with dispersal strength on y-axis and distance (km) on x-axis

x = seq(0, 50, 0.0001)

y = (1/sigma*sqrt(2*3.14))*exp(-(x^2)/2*sigma^2)  #1 is the shape parameter?
y


y = 1/(sigma*sqrt(2))*exp(-(sqrt(2*x))/sigma) #Lapacian

y2 = 2/(sigma*sqrt(2))*exp(-(sqrt(2*x))/sigma) #not sure why this shape is not changing

plot(x, y2, type="l")

plot(x, y, type="l")


par(bg="white")
plot(x,y,type="l", col=NA, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
polygon(x,y, col="grey", border=NA)
