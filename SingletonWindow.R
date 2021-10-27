# LT 14/10/2021

# singletons in a samples are a window on the SFS of the population

# example with annon modal SFS

# population size
N = 1e2

#sample size
n = 10

# back to the basics: how a mutation with AC = i in the population is projected 
# in the smaller sample
i=50
dhyper(0:n, i, N-i, n)

i=1
dhyper(0:n, i, N-i, n)

# now we do the opposite: given a mutation is a singleton in the sample, where
# it come from

i = 1:N
dhyper(rep(1,N), i, N-i, n)
plot(dhyper(rep(1,N), i, N-i, n))
lines(dhyper(rep(1,N), i, N-i, 5))

# Ok generate a figure for a reasonable population size, let's say 1e6
# and range of of samples from 10 to 1000

# we don't need all numbers: think about that
N = 1e4
i = 1:N
n=1000
plot(dhyper(1, i, N-i, N/10), type='l', col=2, xlab= 'Allele count', ylab='Contribution to singletons', log='')
lines(dhyper(rep(1,N), i, N-i, N/100), col=3)
lines(dhyper(rep(1,N), i, N-i, N/1000), col=4)
legend('topright', legend=c("n=1000","n=100","n=10"), col=2:4, lwd=2)


# example of how proportion of singletons changes with sample size and SFS shape

SFS = 1/(1:N)
SFS = SFS / sum(SFS)

modalSFS = dlnorm(1:N, meanlog=2.5, sdlog=0.5)
modalSFS = modalSFS/sum(modalSFS)

lines(SFS)
lines(modalSFS)

# now plot how prop of singletons changes with sample size
ns = c(5000, 1000, 100, 10)
i = 1:N
ps.SFS = sapply(ns, function(n) sum(dhyper(1, i, N-i, n) * SFS)/(1-sum(dhyper(0, i, N-i, n) * SFS)))
ps.modalSFS = sapply(ns, function(n) sum(dhyper(1, i, N-i, n) * modalSFS)/(1-sum(dhyper(0, i, N-i, n) * modalSFS)))
