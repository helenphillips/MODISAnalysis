## Author: Helen Phillips
## Created: 30th July 2014

## Function to create a vector of weighted means the same order as a vector of pixel values


weight_matrix <- function(pixel.tile){
	
	n.row <- sqrt(ncol(pixel.tile))
	n <- floor(n.row/2)


	x <- 1 ## the central pixel, so has a the maxium weight possible
	
	## This for loops creates a vector of the weights we need for each "band" of pixels around the central pixel

	for(i in 1:n){
		weight <- exp(-((i^2)/(n^2)))
		x <- c(x, rep(weight, i*8))
	}

	v <- rev(x)

	## We need the vector in the same order as teh vector of pixels, which at the moment it isn't.
	## If we put it in a matrix we want the first position of the vector to be at the centre and the rest spiraling out
	## this function will do that for us using a recursion
	## HUGE thank you to http://rosettacode.org/wiki/Spiral_matrix#R for putting their code online to do this


	spiralv<-function(v){
		n<-sqrt(length(v))
		if(n!=floor(n)) stop(simpleError("length of v should be a square of an integer"))
		if(n==0) stop(simpleError("v should be of positive length"))
		if(n==1) M<-matrix(v,1,1)
		else M<-rbind(v[1:n],cbind(spiralv(v[(2*n):(n^2)])[(n-1):1,(n-1):1],v[(n+1):(2*n-1)]))
		M
	}


	options("expressions") ## This number needs to be high
	 options(expressions = 10000)


	test <- spiralv(v)
## won't be able to see this in full, head(test) might help though.
# this has created a matrix, where higher numbers indicate that the pixel is further away from teh central pixel
# compared to lower numbers.


	## Then turn it back into a vector, which is now in the correct order

	test2 <- as.vector(matrix(test, nrow = 1))

	return(test2)
}
