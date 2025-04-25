
# Section 6: Incorporating Run Time into Design ---------------------------

# Figure 20

# Functions ---------------------------------------------------------------

# Figure 20 ---------------------------------------------------------------

sig <- function(x,
                ceiling
){
  a <- 5 #5 or 10
  b <- 3
  c <- ceiling
  d <- 0.5
  if (x<=c){
    sigmoidal <- b/(1+exp(-a*(x-c)))
  }
  if (x>c){
    sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
  }
  return(sigmoidal)
}
x <- seq(0,10,by=0.01)
# plot(x,sapply(x,FUN=sig,ceiling=5),
#      type="l",
#      lwd=3,
#      xlab="Total Run Time",
#      ylab="Penalty",
#      )


#png(file="sigmoidal.png",width=2150,height=1315,res=200)
plot(x,sapply(x,FUN=sig,ceiling=5),
     type="l",
     lwd=3,
     xlab="Total Run Time",
     ylab="Penalty",
     xaxt = "n",
     main = "Sigmoidal Function (and a Linear Term)"
)

axis(1, at = c(5),
     labels = c("Ceiling"))

#dev.off()


#larger margins for poster/presentation
#png(file="sigmoidal_larger_margins.png",width=1612.5,height=986.25,res=200)
plot(x,sapply(x,FUN=sig,ceiling=5),
     type="l",
     lwd=3,
     xlab="Total Run Time",
     ylab="Penalty",
     xaxt = "n",
     main = "Sigmoidal Function (and a Linear Term)"
)

axis(1, at = c(5),
     labels = c("Ceiling"))

#dev.off()


#larger margins for poster/presentation
#changed x axis label
#png(file="new_sigmoidal_larger_margins.png",width=1612.5,height=986.25,res=200)
plot(x,sapply(x,FUN=sig,ceiling=5),
     type="l",
     lwd=3,
     xlab="Expected Total Run Time (Tc)",
     ylab="Penalty",
     xaxt = "n",
     main = "Sigmoidal Function (and a Linear Term)"
)

axis(1, at = c(5),
     labels = c("Ceiling"))

#dev.off()


# End of Section 6 --------------------------------------------------------
