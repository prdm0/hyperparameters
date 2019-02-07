S2 = function(y1, y2){
   D = outer(y1, y2, '-')
   D = D^2
   D_no_zero = D[which(!D == 0)]
   median(D_no_zero)
}
# Example: 

set.seed(0)
size_vector <- 2
y1 <- runif(n = size_vector, min = 0, max = 1)
y2 <- runif(n = size_vector, min = 0, max = 1)

S2(y1 = y1, y2 = y2)
S2_cpp(y1 = y1, y2 = y2)