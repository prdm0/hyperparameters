S2 = function(y1, y2){
   D = outer(y1, y2, '-')
   D = D^2
   D_no_zero = D[which(!D == 0)]
   median(D_no_zero)
}