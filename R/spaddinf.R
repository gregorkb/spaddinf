# functions for sparse additive inference:

#' Plot an image of a matrix
implot <- function(x,legend=TRUE)
{
	x.new <- t(x)[,dim(x)[1]:1]
	n <- 100

	if(min(x.new)*max(x.new)<0)
	{
		into.red <- round(n*min(abs(min(x.new))/max(x.new),1)) # deepest red
		into.blue <- round(n*min(max(x.new)/abs(min(x.new)),1))

		hsvvector <- c(hsv(s=into.red:0/n,1,1),hsv(s=1:into.blue/n,2/3,1))

	}else hsvvector <- hsv(s=1:n/n,2/3,1)

	image(x.new,col = hsvvector ,xaxt="n" , yaxt="n",bty="n")
	box()

	# Ersatz legend:
	if(legend==TRUE)
	{
		hsvvector.legend <- c(hsv(s=50:0/50,1,1),hsv(s=1:50/50,2/3,1))
		x.pos <- .6
		y.pos <- .9
		points(seq(x.pos,x.pos+.3,length=101),rep(y.pos,101),col=hsvvector.legend ,cex=2.5,pch=15)
		text(x.pos,y.pos+.01,"-",pos=3,cex=1)
		text(x.pos+.3,y.pos+.01,"+",pos=3,cex=1)
	}
	# zero should be white, but negative numbers should be red...
}


#################### Functions for building design matrices of spline evaluations
make.HQ <- function (X, d, degree,knots.list=list())
{

	if( length(knots.list) == 0 ){

		no.input.knots <- 1

	} else { no.input.knots <- 0}

    q <- ncol(X)
    n <- nrow(X)
    H <- H.tilde <- matrix(0, n, q * d)
    Q.inv <- matrix(0, q * d, q * d)
    for (j in 1:q) {
    	if(no.input.knots)
    	{
        	knots <- quantile(X[, j], seq(0, 1, length = d - degree + 1))

        	knots.list[[length(knots.list)+1]] <- knots

        } else {

        	knots <- knots.list[[j]]
        }

        boundary.knots <- range(knots)
        new.knots <- sort(c(rep(boundary.knots, degree), knots))
        x <- X[, j]
        bspline.mat <- spline.des(x, knots = new.knots, outer.ok = TRUE,derivs = rep(0, length(x)), ord = degree + 1)$design
        H[, ((j - 1) * d + 1):(d * j)] <- bspline.mat
        Q.s <- chol(t(bspline.mat) %*% bspline.mat/n)
        Q.inv[((j - 1) * d + 1):(d * j), ((j - 1) * d + 1):(d * j)] <- solve(Q.s)
        H.tilde[, ((j - 1) * d + 1):(d * j)] <- bspline.mat %*% solve(Q.s)
    }

    output <- list(H = H, Q.inv = Q.inv, H.tilde = H.tilde, knots.list = knots.list)
    return(output)
}


make.H <- function (X, d, degree, knots.list = list())
{
    if (length(knots.list) == 0) {

        no.input.knots <- 1
    }
    else{

        no.input.knots <- 0
    }

    q <- ncol(X)
    n <- nrow(X)
    H <- matrix(0, n, q * d)

    for(j in 1:q)
    {
        if(no.input.knots) {
            knots <- quantile(X[, j], seq(0, 1, length = d - degree + 1))
            knots.list[[length(knots.list) + 1]] <- knots
        }
        else {

            knots <- knots.list[[j]]
        }

        boundary.knots <- range(knots)
        new.knots <- sort(c(rep(boundary.knots, degree), knots))
        x <- X[, j]
        bspline.mat <- spline.des(x, knots = new.knots, outer.ok = TRUE,derivs = rep(0, length(x)), ord = degree + 1)$design
        H[, ((j - 1) * d + 1):(d * j)] <- bspline.mat

    }
    output <- list(H = H, knots.list = knots.list)
    return(output)
}


make.W <- function(X1,d.re,degree.re)
{

	knots <- quantile(X1,seq(0,1,length= d.re-degree.re +1))
	# knots <- seq(-2.5,2.5,length= d.re-degree.re +1)
	boundary.knots <- range(knots)
	new.knots <- sort(c(rep(boundary.knots, degree.re),knots))

	W <- spline.des(X1, knots = new.knots,outer.ok=TRUE,derivs=rep(0,length(X1)),ord=degree.re +1)$design

	W.cent <- scale(W,center=TRUE,scale=FALSE)

	return(W)

}

getX <- function(n,q,r)
{

	R <- r^abs( outer(1:q,1:q,"-"))
	P <- 2*sin( R * pi / 6)

	X <- (pnorm( matrix(rnorm(n* q),ncol= q) %*% chol(P)) - .5) * 5

	output <- list(	X = X,
					        r = r)

	return(output)

}

#' Helps build the basis functions for Legendre polynomials
pcws.poly <- function(x,left.endpoints,K){

  vec <- numeric(length(left.endpoints)*(K+1))
  int <- sum(x >= left.endpoints) # in which interval does it fall?
  vec[((int-1)*(K+1)+1):((int)*(K+1))] <- 1*(x - left.endpoints[int])^c(0:K)
  return(vec)

}

#' Fit the desparsified lasso presmoothing estimator with cubic B-splines
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param lambda the tuning parameter for fitting the group lasso estimate for the bias correction
#' @param eta the tuning parameter for the group lasso projection of one set of basis functions onto those of the other covariates.
#' @param n.foi the number of functions (first columns of \code{X} on which to compute the desparsified lasso presmoothing estimator.
#' @returns a list with the fitted functions etc.
#'
#' @examples
#' # set parameters
#' r <- .5
#' n <- 100
#' q <- 10
#' d.pre <- max(30,floor(2 * sqrt(n)))
#' n.foi <- 6
#' lambda <- 4
#' eta <- 3
#'
#' f1 <- function(x){-sin(x*2)}
#' f2 <- function(x){x^2 - 25/12}
#' f3 <- function(x){x}
#' f4 <- function(x){exp(-x)-2/5*sinh(5/2)}
#' f5 <- function(x){x*0}
#' f6 <- function(x){x*0}
#'
#' # generate data
#' X <- getX(n,q,r)$X
#' signal.uncent <- cbind( f1(X[,1]),f2(X[,2]),f3(X[,3]),f4(X[,4]),f5(X[,5]),f6(X[,6]) )
#' signal.means <- apply(signal.uncent,2,mean)
#' signal <- signal.uncent - matrix(signal.means,n,6,byrow=TRUE)
#' noise <- rnorm(n)
#' Y <- apply(signal,1,sum) + noise - mean(noise)
#' spadd.presmth.Bspl.out <- spadd.presmth.Bspl(X,Y,d.pre,lambda,eta,n.foi)
#' f.hat <- spadd.presmth.Bspl.out$f.hat
#'
#' # plot output
#' par(mfrow=c(2,3))
#' for( j in 1:6){
#'
#'   plot(f.hat[[j]],min(X[,j]),max(X[,j]))
#'
#'   x.seq <- seq(min(X[,j]),max(X[,j]),length=300)
#'   fj <- f[[j]](x.seq) - mean(f[[j]](X[,j])) # must center the true function in order to make comparisons
#'   lines(fj~x.seq,lty=2)
#'
#' }
spadd.presmth.Bspl <- function(X,Y,d.pre,lambda,eta,n.foi)
{

  q <- ncol(X)
  n <- nrow(X)

  # make cubic B-splines basis function design
  HH <- HH.tilde <- matrix(NA,n,0)
  groups <- numeric()
  QQ.inv <- vector("list",length=q)
  knots.list <- vector("list",length=q)
  emp.cent <- vector("list",length=q)

  for( j in 1:q )
  {

    int.knots <- quantile(X[,j],seq(0,1,length=d.pre-2+1)) # add one, so that one can be removed after centering to restore full-rank.
    boundary.knots <- range(int.knots)
    all.knots <- sort(c(rep(boundary.knots,3),int.knots))
    knots.list[[j]] <- all.knots

    Bj <- spline.des(all.knots,X[,j],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
    emp.cent[[j]] <- apply(Bj,2,mean)
    Bj.cent <- Bj - matrix(emp.cent[[j]],n,d.pre,byrow=TRUE)

    # construct matrix in which l2 norm of function is a quadratic form
    M <- t(Bj.cent) %*% Bj.cent / n

    Q <- chol(M)
    Q.inv <- solve(Q)
    QQ.inv[[j]] <- Q.inv

    # construct basis function matrices
    HH.tilde <- cbind(HH.tilde,Bj.cent %*% Q.inv)
    HH <- cbind(HH,Bj.cent)
    groups <- c(groups,rep(j,d.pre))

  }

  # get the group lasso estimators
  grplasso.out <- grplasso(	y = Y,
                            x = HH.tilde,
                            index = groups,
                            lambda = lambda,
                            model = LinReg(),
                            center = FALSE,
                            standardize = FALSE,
                            control = grpl.control(trace=0))

  beta.tilde <- grplasso.out$coef
  f.lasso <- grplasso.out$fitted
  selected <- grplasso.out$norms.pen != 0

  f.lasso.design <- f.hat.design <- matrix(NA,nrow(X),n.foi)
  f.hat <- list()

  for(j in 1:n.foi)
  {

    # Now get Z1, the matrix replacing the projection of H_1 onto the
    # orthogonal complement of the other columns of H.  Use the group
    # lasso to produce a projection having orthogonality
    # controlled by lambda.
    ind <- which(groups == j)

    Zj <- matrix(0,n,d.pre)

    for(l in 1:d.pre)
    {
      Z.model <- grplasso(y = HH[, ind[l]],
                          x = HH.tilde[,-ind],
                          index = groups[-ind],
                          lambda = eta,
                          model = LinReg(),
                          center = FALSE,
                          standardize = FALSE,
                          control = grpl.control( trace = 0 )
      )

      Zj[,l] <- HH[,ind[l]] - Z.model$fitted

    }

    # Construct the presmoothing estimator f.1.hat of f_1 based on the
    # nonparametric desparsified Lasso

    Xj <- HH[,ind]
    fj.lasso <- HH.tilde[,ind] %*% beta.tilde[ind] # grp lasso estimate of f.1
    Yj.lasso <- Y - HH.tilde[,-ind] %*% beta.tilde[-ind]

    betaj.hat <- solve( t(Zj) %*% Xj ) %*% t(Zj) %*% Yj.lasso
    fj.hat <- Xj %*% betaj.hat	# first-stage estimator of f.1

    f.lasso.design[,j] <- fj.lasso
    f.hat.design[,j] <- fj.hat

    # export actual function estimate
    f.hat[[j]] <- eval(parse(text=paste("function(x)
                                        {

                                        x <- round(x,10)
                                        x.mat <- spline.des(",paste("c(",paste(round(knots.list[[j]],6),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
                                        x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent[[j]],collapse=","),"),length(x),",d.pre,sep=""),",byrow=TRUE)
                                        f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(betaj.hat,collapse=","),")",sep=""),")
                                        return(f.hat)

                                        }"
      			)))
  }

  output <- list(	f.lasso = f.lasso,
                  f.lasso.design = f.lasso.design,
                  f.hat.design = f.hat.design,
                  f.hat = f.hat)

	return(output)

}

#' Fit the desparsified lasso presmoothing estimator with Legendre polynomials
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param lambda the tuning parameter for fitting the group lasso estimate for the bias correction
#' @param eta the tuning parameter for the group lasso projection of one set of basis functions onto those of the other covariates.
#' @param n.foi the number of functions (first columns of \code{X} on which to compute the desparsified lasso presmoothing estimator.
#' @param K the order of the Legendre polynomials. E.g. \code{K=0} fits piecwise constant, \code{K=1} fits piecewise linear functions.
#' @returns a list with the fitted functions etc.
#'
#' @examples
#' r <- .5
# n <- 100
# q <- 20
# d.pre <- max(20,floor(2 * sqrt(n)))
# n.foi <- 6
# lambda <- .5
# eta <- 3
#
# f <- vector("list",q)
# f[[1]] <- function(x){-sin(x*2)}
# f[[2]] <- function(x){x^2 - 25/12}
# f[[3]] <- function(x){x}
# f[[4]] <- function(x){exp(-x)-2/5*sinh(5/2)}
# f[[5]] <- function(x){x*0}
# f[[6]] <- function(x){x*0}
#
# ######## Generate data
# X <- getX(n,q,r)$X
# signal.uncent <- cbind( f[[1]](X[,1]),f[[2]](X[,2]),f[[3]](X[,3]),f[[4]](X[,4]),f[[5]](X[,5]),f[[6]](X[,6]) )
# signal.means <- apply(signal.uncent,2,mean)
# signal <- signal.uncent - matrix(signal.means,n,6,byrow=TRUE)
# noise <- rnorm(n)
# Y <- apply(signal,1,sum) + noise - mean(noise)
#
# spadd.presmth.Legr.out <- spadd.presmth.Legr(X,Y,d.pre=10,lambda=.1,eta,n.foi,K=1)
#
#
# f.hat <- spadd.presmth.Legr.out$f.hat
# f.hat.design <- spadd.presmth.Legr.out$f.hat.design
# left.endpoints.list <- spadd.presmth.Legr.out$left.endpoints.list
# ######## plot output
#
# par(mfrow = c(2,3),
#     oma = c(2.1,3,2.1,.1),
#     mar = c(0,0,0,0),
#     lwd = 1,
#     cex.axis = .8,
#     cex.lab = .85,
#     tck = -.01)
#
# for( j in 1:6)
# {
#
#   plot(NA,
#        xlim=range(X[,j]),
#        ylim=range(f.hat.design[,j],f[[j]](X[,j])),
#        xlab="",
#        ylab="",
#        xaxt="n",
#        yaxt="n")
#
#   # points(f.hat.design[order(X[,j]),j]~sort(X[,j]))
#   # plot each fitted Legendre polynomial function
#   for(l in 1:(length(left.endpoints.list[[j]])-1))
#   {
#
#     x1 <- left.endpoints.list[[j]][l]
#     x2 <- left.endpoints.list[[j]][l+1]
#     plot(f.hat[[j]],x1,x2-1e-4,col=rgb(.545,0,0,1),lwd=1.5,add=TRUE)
#   }
#
#   plot(f.hat[[j]],x2,max(X[,j])-1e-4,col=rgb(.545,0,0,1),lwd=1.5,add=TRUE)
#
#   x.seq <- seq(min(X[,j]),max(X[,j]),length=300)
#   fj <- f[[j]](x.seq) - mean(f[[j]](X[,j])) # must center the true function in order to make comparisons
#   lines(fj~x.seq,lty=2)
#
#   if(j >=4) axis(1,padj=-1.5)
#   if(j ==1 | j==4) axis(2,padj=1)
#
#   x.pos <- grconvertX(.5,from="nfc",to="user")
#   y.pos <- grconvertY(.08,from="nfc",to="user")
#
#   # text(x = x.pos, y = y.pos, labels = paste("coverage of f_",j,sep=""))
#   if(j == 1) text(x = x.pos, y = y.pos, labels = bquote(f[1](x)== -sin(2*x) ))
#   if(j == 2) text(x = x.pos, y = y.pos, labels = bquote(f[2](x)== x^2 - 25/12 ))
#   if(j == 3) text(x = x.pos, y = y.pos, labels = bquote(f[3](x)== x ))
#   if(j == 4) text(x = x.pos, y = y.pos, labels = bquote(f[4](x)== exp(-x) - (2/5) * sinh(5/2) ))
#   if(j == 5) text(x = x.pos, y = y.pos, labels = bquote(f[5](x)== 0))
#   if(j == 6) text(x = x.pos, y = y.pos, labels = bquote(f[6](x)== 0 ))
#
#   if(j==6) mtext(side=1,bquote(x),line=1.5,cex=.8)
#
# }
#
# x.pos <- grconvertX(.075,from="ndc",to="user")
# y.pos <- grconvertY(.1,from="ndc",to="user")
#
# mtext(side=3,outer=TRUE,text=paste("Presmoother with Legendre Polynomials, n = ",n,", d = ", d.pre,sep=""),line=.5,cex=.8)
spadd.presmth.Legr <- function(X,Y,d.pre,lambda,eta,n.foi,K=1)
{
  q <- ncol(X)
  n <- nrow(X)

  # make cubic B-splines basis function design
  HH <- HH.tilde <- matrix(NA,n,0)
  groups <- numeric()
  QQ.inv <- vector("list",length=q)
  left.endpoints.list <- vector("list",length=q)
  emp.cent <- vector("list",length=q)

  for( j in 1:q )
  {

    left.endpoints <- quantile(X[,j],seq(0,1,length=d.pre+2))[-(d.pre+2)] # add one interval so that one can be removed after centering to restore full-rank.
    left.endpoints.list[[j]] <- left.endpoints

    # now make the Legendre polynomial basis
    Bj <- t(sapply(X[,j],FUN=pcws.poly,left.endpoints = left.endpoints,K=K))[,-c(1:(K+1))] # remove columns corresponding to first interval
    emp.cent[[j]] <- apply(Bj,2,mean)
    Bj.cent <- Bj - matrix(emp.cent[[j]],n,d.pre*(K+1),byrow=TRUE)

    # construct matrix in which l2 norm of function is a quadratic form
    M <- t(Bj.cent) %*% Bj.cent / n

    Q <- chol(M)
    Q.inv <- solve(Q)
    QQ.inv[[j]] <- Q.inv

    # construct basis function matrices
    HH.tilde <- cbind(HH.tilde,Bj.cent %*% Q.inv)
    HH <- cbind(HH,Bj.cent)
    groups <- c(groups,rep(j,d.pre*(K+1)))

  }

  # get the group lasso estimators
  grplasso.out <- grplasso(	y = Y,
                            x = HH.tilde,
                            index = groups,
                            lambda = lambda,
                            model = LinReg(),
                            center = FALSE,
                            standardize = FALSE,
                            control = grpl.control(trace=0))

  beta.tilde <- grplasso.out$coef
  f.lasso <- grplasso.out$fitted
  selected <- grplasso.out$norms.pen != 0

  f.lasso.design <- f.hat.design <- matrix(NA,nrow(X),n.foi)
  f.hat <- list()

  for(j in 1:n.foi)
  {

    # Now get Z1, the matrix replacing the projection of H_1 onto the
    # orthogonal complement of the other columns of H.  Use the group
    # lasso to produce a projection having orthogonality
    # controlled by lambda.
    ind <- which(groups == j)

    Zj <- matrix(0,n,d.pre*(K+1))

    for(l in 1:(d.pre*(K+1)))
    {
      Z.model <- grplasso(y = HH[, ind[l]],
                          x = HH.tilde[,-ind],
                          index = groups[-ind],
                          lambda = eta,
                          model = LinReg(),
                          center = FALSE,
                          standardize = FALSE,
                          control = grpl.control( trace = 0 )
      )

      Zj[,l] <- HH[,ind[l]] - Z.model$fitted

    }

    # Construct the presmoothing estimator f.1.hat of f_1 based on the
    # nonparametric desparsified Lasso
    Xj <- HH[,ind]
    fj.lasso <- HH.tilde[,ind] %*% beta.tilde[ind] # grp lasso estimate of f.1
    Yj.lasso <- Y - HH.tilde[,-ind] %*% beta.tilde[-ind]

    betaj.hat <- solve( t(Zj) %*% Xj ) %*% t(Zj) %*% Yj.lasso
    fj.hat <- Xj %*% betaj.hat	# first-stage estimator of f.1

    f.lasso.design[,j] <- fj.lasso
    f.hat.design[,j] <- fj.hat

    # export actual function estimate
    f.hat[[j]] <- eval(parse(text=paste("function(x)
                                        {
                                        x <- round(x,10)
                                        x.mat <- t(sapply(x,FUN=pcws.poly,left.endpoints = c(",paste(round(left.endpoints.list[[j]],10),collapse=","),"), K=",K,"))[,-c(1:(",K,"+1))]
                                        x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent[[j]],collapse=","),"),length(x),",d.pre*(K+1),sep=""),",byrow=TRUE)
                                        f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(betaj.hat,collapse=","),")",sep=""),")
                                        return(f.hat)
                                        }"
        			)))
  }

  output <- list(	f.lasso = f.lasso,
                  f.lasso.design = f.lasso.design,
                  f.hat.design = f.hat.design,
                  f.hat = f.hat,
                  left.endpoints.list = left.endpoints.list)

  return(output)

}

#' Fit simple nonparametric regression model with cubic B-splines
#'
#' @param Y a response vector
#' @param X vector of covariate observations
#' @param d the number functions in the cubic B-spline basis
#' @return a list containing the fitted function and a vector containing the values of the fitted function at the design points
#'
#' @examples
#' # generate some data
#' n <- 500
#' d <- 10
#' f <- function(x){-sin(x*2)}
#' X <- runif(n,-2.5,2.5)
#' noise <- rnorm(n)
#' Y <- f(X) - mean(f(X)) + noise - mean(noise)
#'
#' # get nonparametric estimate of f
#' smth.Bspl.out <- smth.Bspl(Y,X,d)
#' f.hat <- smth.Bspl.out$f.hat
#'
#' # plot results
#' par(mfrow=c(1,1))
#' plot(Y~X)
#' plot(f.hat,min(X)+1e-3,max(X)-1e-3,add=TRUE,lwd=1.5,col=rgb(0,0,.545))
#' x.seq <- seq(min(X),max(X),length=500)
#' f.cent <- f(x.seq) - mean(f(X))
#' lines(f.cent~x.seq,lty=2)
smth.Bspl <- function(Y,X,d)
{

  int.knots <- quantile(X,seq(0,1,length=d-2+1)) # add one, so that one can be removed after centering to restore full-rank.
  boundary.knots <- range(int.knots)
  all.knots <- sort(c(rep(boundary.knots,3),int.knots))

  B <- spline.des(all.knots,X,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
  emp.cent <- apply(B,2,mean)
  B.cent <- B - matrix(emp.cent,n,d,byrow=TRUE)

  beta.hat <- as.numeric(solve( t(B.cent) %*% B.cent) %*% t(B.cent) %*% (Y - mean(Y)))
  f.hat.design <- as.numeric(B.cent %*% beta.hat)

  f.hat <- eval(parse(text=paste("function(x){
                                 x <- round(x,10)
                                 x.mat <- spline.des(",paste("c(",paste(round(all.knots,6),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
                                 x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent,collapse=","),"),length(x),",d,sep=""),",byrow=TRUE)
                                 f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(beta.hat,collapse=","),")",sep=""),")
                                 return(f.hat)
}"
        			)))

  output <- list(f.hat.design = f.hat.design,
                 f.hat = f.hat,
                 d = d)

  return(output)

}


#################### crossvalidation functions
cv.B.splines  <- function (X, Y, d.seq, K, degree, std.err = 1, plot = FALSE)
{
    n <- length(X)
    cv.msep <- cv.msep.se <- numeric()
    for (l in 1:length(d.seq)) {
        pred <- numeric(n)
        fold.msep <- numeric(K)
        for (k in 1:K) {
            F.k <- k + (1:(floor(n/K) + 1) - 1) * K
            F.k <- F.k[F.k <= n]
            knots <- quantile(X[-F.k], seq(0, 1, length = d.seq[l] -
                degree + 1))
            boundary.knots <- range(knots)
            new.knots <- sort(c(rep(boundary.knots, degree),
                knots))
            W.k.train <- spline.des(X[-F.k], knots = new.knots,
                outer.ok = TRUE, derivs = rep(0, length(X[-F.k])),
                ord = degree + 1)$design
            W.k.test <- spline.des(X[F.k], knots = new.knots,
                outer.ok = TRUE, derivs = rep(0, length(X[F.k])),
                ord = degree + 1)$design
            pred[F.k] <- W.k.test %*% solve(t(W.k.train) %*%
                W.k.train) %*% t(W.k.train) %*% Y[-F.k]
            fold.msep[k] <- mean((pred[F.k] - Y[F.k])^2)
        }
        cv.msep[l] <- mean((pred - Y)^2)
        cv.msep.se[l] <- sd(fold.msep)/sqrt(K)
    }
    min.cv.msep <- cv.msep[which(cv.msep == min(cv.msep))]
    cv.msep.se.at.min <- cv.msep.se[which(cv.msep == min(cv.msep))]
    d.cv.std.err <- d.seq[max(which(cv.msep <= min.cv.msep + std.err * cv.msep.se.at.min))]
    d.cv.min <- d.seq[max(which(cv.msep == min.cv.msep))]
    output <- list(cv.msep = cv.msep, cv.msep.se = cv.msep.se,
        d.cv.std.err = d.cv.std.err, d.cv.min = d.cv.min)
    if (plot == TRUE) {
        plot(cv.msep ~ d.seq, ylim = range(c(cv.msep, cv.msep + std.err * cv.msep.se)))
        lines(cv.msep + cv.msep.se ~ d.seq)
    }
    return(output)
}

cv.loc.pol <- function(X,Y,K,Ker,r,num.h,plot=FALSE)
{

	n <- length(X)
	h.min <- max(diff(sort(X)))/3
	h.max <- diff(range(X))/4
	h.seq <- seq(h.min^(1/3),h.max^(1/3),length=num.h)^3

	cv.msep <- cv.msep.se <- numeric()

	for( l in 1:length(h.seq))
	{

		pred <- numeric(n)
		fold.msep <- numeric(K)

		for(k in 1:K)
		{

			F.k <- k + (1:(floor(n/K)+1) - 1)*K
		    F.k <- F.k[F.k<=n]

			data.no.F.k <- list(X = X[-F.k],Y=Y[-F.k])
		    pred[F.k] <- sapply(X[F.k],FUN=loc.pol.x,data=data.no.F.k,h=h.seq[l],Ker=Ker,r=r)
			fold.msep[k] <- mean( (pred[F.k] - Y[F.k])^2 )


		}

		cv.msep[l] <- mean( (pred - Y)^2)
		cv.msep.se[l] <- sd(fold.msep)/sqrt(K)

  	}

	min.cv.msep <- cv.msep[which(cv.msep == min(cv.msep))]
	cv.msep.se.at.min <- cv.msep.se[which(cv.msep == min(cv.msep))]
	h.cv <- h.seq[min(which(cv.msep <= min.cv.msep + 1*cv.msep.se.at.min))]

	output <- list(	cv.msep = cv.msep,
					cv.msep.se = cv.msep.se,
					h.cv = h.cv)


	if(plot==TRUE)
	{

		plot(cv.msep~h.seq,ylim=range(c(cv.msep, cv.msep + 1*cv.msep.se)))
		lines(cv.msep + cv.msep.se~h.seq)

	}

	return(output)

}


################################################################################################################
#################### Selection of lambda and eta
################################################################################################################

spaddinf.presmth.cv <- function(X,Y,fois,d.pre,degree.pre,length.lambda.seq,length.eta.seq,K.lambda,K.eta,std.errs = .05,plot=FALSE)
{

		grps <- c(1:ncol(X)) %x% rep(1,d.pre)

		##############################################################################
		##############################################################################
		###### Evaluate spline basis functions at design points
		###### and store in new design matrix H for j = 1,...,q
		##############################################################################
		##############################################################################

		HQ <- make.HQ(X,d.pre,degree.pre)
		H <- HQ$H
		H.tilde <- HQ$H.tilde
		Q.inv <- HQ$Q.inv

		lambda.max <- lambdamax(y = Y,
								x = H.tilde,
								index = grps,
								model = LinReg(),
						        center = FALSE,
						        standardize = FALSE,
						        control = grpl.control(trace=0))

		lambda.exp <- exp(seq(log(lambda.max),0,length = length.lambda.seq)) - 1
		lambda.seq <- lambda.exp * lambda.max/max(lambda.exp)

		indscramble <- sample(1:n,n)

		h <- 1

		MSEP.fold.lambda <- matrix(0,K.lambda,length.lambda.seq)

		for(k in 1:K.lambda)
		{

			fold.ind <- indscramble[h:(h + (floor(n/K.lambda)+1)*(k<=(n%%K.lambda)) + floor(n/K.lambda)*(k>(n%%K.lambda)) - 1 )]
			h <- h + (floor(n/K.lambda)+1)*(k<=(n%%K.lambda)) + floor(n/K.lambda)*(k>(n%%K.lambda))
			n.fold <- length(fold.ind)

			X.fold.train <- X[-fold.ind,]
			X.fold.test <- X[fold.ind,]
			Y.fold.train <- Y[-fold.ind]
			Y.fold.test <- Y[fold.ind]

			HQ.fold.train <- make.HQ(X.fold.train,d.pre,degree.pre)
			H.fold.train <- HQ.fold.train$H
			H.tilde.fold.train <- HQ.fold.train$H.tilde
			Q.inv.fold.train <- HQ.fold.train$Q.inv

			H.fold.test <- make.H(X.fold.test,d.pre,degree.pre,knots.list=HQ.fold.train$knots.list)$H

			grplasso.out.fold <- grplasso(	y = Y.fold.train,
											x = H.tilde.fold.train,
											index = grps,
											# lambda = lambda.seq,
											lambda = lambda.seq*(K.lambda-1)/K.lambda,
											model = LinReg(),
						                    center = FALSE,
						                    standardize = FALSE,
						                    control = grpl.control(trace=0))

			beta.lasso.fold <- Q.inv.fold.train %*% grplasso.out.fold$coef
			f.lasso.fold <- H.fold.test %*% beta.lasso.fold

			MSEP.fold.lambda[k,] <- apply( (f.lasso.fold - Y.fold.test %*% matrix(1,1, length.lambda.seq))^2,2,mean)

			print(paste("lambda fold: ", k ,sep=""))

		}

		cv.MSEPs.lambda <- apply(MSEP.fold.lambda,2,mean)
		cv.MSEP.se.lambda  <- apply(MSEP.fold.lambda,2,sd)*(sqrt((K.lambda-1)/K.lambda))/sqrt(K.lambda)

		# use one standard error rule from Tibshirani's notes (or rather std.errs standard errors)
		# err on the side of smaller or larger lambda? greater methinks.

		min.cv.MSEPs.lambda <- cv.MSEPs.lambda[which(cv.MSEPs.lambda==min(cv.MSEPs.lambda))]
		min.cv.MSEP.se.lambda <- cv.MSEP.se.lambda[which(cv.MSEPs.lambda==min(cv.MSEPs.lambda))]
		which.lambda <- min(which(cv.MSEPs.lambda <= min.cv.MSEPs.lambda + std.errs * min.cv.MSEP.se.lambda))

		cv.lambda <- lambda.seq[which.lambda]

	##############################################################################
	##############################################################################
	## END-CROSSVALIDATION for lambda
	##############################################################################
	##############################################################################

	##############################################################################
	##############################################################################
	## CROSSVALIDATION for choosing eta for the Lasso procedure producing Z.j
	##############################################################################
	##############################################################################

	# Do only for j = 1

    j <- 1

    H.ind.j <- ((fois[j]-1)*d.pre + 1):(fois[j]*d.pre)

	eta.max <- numeric()

	for(l in 1:d.pre)
	{
		eta.max[l] <- lambdamax(	y = H[, H.ind.j[l]],
									x = H.tilde[,-H.ind.j],  # use H.tilde to impose the right penalty!
									index = grps[-H.ind.j],
									model = LinReg(),
				                    center = FALSE,
				                    standardize = FALSE,
				                    control = grpl.control(trace=0))
	}

	eta.exp <- exp(seq(log(max(eta.max,1.05)), 0, length = length.eta.seq)) - 1
	eta.seq <- eta.exp * max(eta.max) / max(eta.exp)

	indscramble <- sample(1:n,n)

	PRED.eta <- array(0,dim=c(n,d.pre,length.eta.seq))

	h <- 1

	MSEP.fold.eta <- matrix(0,K.eta,length.eta.seq)

	for(k in 1:K.eta)
	{

		print(paste("eta fold: ", k ,sep=""))

		fold.ind <- indscramble[h:(h + (floor(n/K.eta)+1)*(k<=(n%%K.eta)) + floor(n/K.eta)*(k>(n%%K.eta)) - 1 )]
		h <- h + (floor(n/K.eta)+1)*(k<=(n%%K.eta)) + floor(n/K.eta)*(k>(n%%K.eta))
		n.fold <- length(fold.ind)

		X.fold.train <- X[-fold.ind,]
		X.fold.test <- X[fold.ind,]

		HQ.fold.train <- make.HQ(X.fold.train,d.pre,degree.pre)
		H.fold.train <- HQ.fold.train$H
		H.tilde.fold.train <- HQ.fold.train$H.tilde
		Q.inv.fold.train <- HQ.fold.train$Q.inv

		H.fold.test <- make.H(X.fold.test,d.pre,degree.pre,knots.list=HQ.fold.train$knots.list)$H

		for(l in 1:d.pre)
		{

			grplasso.out.fold <- grplasso(	y = H.fold.train[,H.ind.j[l]],
											x = H.tilde.fold.train[,-H.ind.j],
											index = grps[-H.ind.j],
											lambda = eta.seq*(K.eta-1)/K.eta,
											model = LinReg(),
							                center = FALSE,
							                standardize = FALSE,
							                control = grpl.control(trace=0))

			gamma.lasso.fold.l <- Q.inv.fold.train[-H.ind.j,-H.ind.j] %*% grplasso.out.fold$coef

			H.lasso.fold.l <- H.fold.test[,-H.ind.j] %*% gamma.lasso.fold.l

			PRED.eta[fold.ind,l,] <- H.lasso.fold.l



		}

		MSEP.fold.eta[k,] <- apply(  (PRED.eta[fold.ind,,] - array(H.fold.test[,H.ind.j],dim=c(n.fold,d.pre,length.eta.seq)))^2 , c(3) , sum)	/(n.fold*d.pre)


	}


	cv.MSEPs.eta <- apply(MSEP.fold.eta,2,mean)
	cv.MSEP.se.eta <- apply(MSEP.fold.eta,2,sd)/sqrt(K.eta)

	min.cv.MSEPs.eta <- cv.MSEPs.eta[which(cv.MSEPs.eta==min(cv.MSEPs.eta))]
	min.cv.MSEP.se.eta <- cv.MSEP.se.eta[which(cv.MSEPs.eta==min(cv.MSEPs.eta))]
	which.eta <- min(which(cv.MSEPs.eta <= min.cv.MSEPs.eta + std.errs* min.cv.MSEP.se.eta))

	cv.eta <- eta.seq[which.eta]



	if(plot==TRUE)
	{

		par(mfrow=c(1,2))

		plot(MSEP.fold.lambda[1,]~lambda.seq,ylim=range(MSEP.fold.lambda))
		for(k in 2:K.lambda) points(MSEP.fold.lambda[k,]~lambda.seq)

		points(cv.MSEPs.lambda~lambda.seq,ylim=range(cv.MSEPs.lambda, cv.MSEPs.lambda + cv.MSEP.se.lambda),pch=19)
		lines(cv.MSEPs.lambda + cv.MSEP.se.lambda ~ lambda.seq)
		abline(v = cv.lambda)


		plot(MSEP.fold.eta[1,]~eta.seq,ylim=range(MSEP.fold.eta))
		for(k in 2:K.eta) points(MSEP.fold.eta[k,]~eta.seq)

		points(cv.MSEPs.eta~eta.seq,pch=19)
		lines(cv.MSEPs.eta + cv.MSEP.se.eta~eta.seq)

		abline(v = cv.eta)

	}

	output <- list( cv.lambda = cv.lambda,
					cv.eta = cv.eta)

	return(output)

}



################################################################################################################
#################### B spline smoother
################################################################################################################

B.spl <- function(x,X,Y,d,degree,A=NULL,boot.band = FALSE,e.for.boot.band = NULL,alpha=0.05)
{

	#A is like the projection matrix for a presmoothing estimator
	if(length(A)==0){A <- diag(length(X))}

	n <- length(X)

	knots <- quantile(X,seq(0,1,length=d-degree+1))
	boundary.knots <- range(knots)
	new.knots <- sort(c(rep(boundary.knots,degree),knots))

	W.x <- spline.des(x = x, knots = new.knots,outer.ok=TRUE,derivs=rep(0,length(x)),ord=degree+1)$design
	W <- spline.des(x = X, knots = new.knots,outer.ok=TRUE,derivs=rep(0,length(X)),ord=degree+1)$design

	beta <- solve(t(W) %*% W) %*% t(W) %*%  Y
	f.B.spl <- W %*% beta
	f.B.spl.x <- W.x %*% beta

	f.B.spl.x.se <- sqrt(diag(W.x %*% solve(t(W) %*% W) %*% t(W) %*% A %*% t(A) %*% W %*% solve(t(W) %*% W) %*% t(W.x)))

	t.star <- NULL

	if(boot.band == TRUE)
	{
		if(length(e.for.boot.band)==0)
		{

			e.for.boot.band <- Y - f.B.spl

		}

		B <- 10000

		wboot.e <- matrix( rnorm(n*B,0,sqrt(rep(e.for.boot.band^2,B))),nrow=n)

		delta.wboot.x <- diag(1/f.B.spl.x.se) %*% W.x %*% solve(t(W) %*% W) %*% t(W) %*% A %*%  wboot.e

		delta.wboot.sup.over.x <- apply(abs(delta.wboot.x),2,max)

		t.star <- quantile(abs(delta.wboot.sup.over.x),1-alpha)

	}

	output <- list(	f.B.spl = f.B.spl,
					f.B.spl.x = f.B.spl.x,
					f.B.spl.x.se = f.B.spl.x.se,
					W = W,
					x = x,
					t.star = t.star,
					alpha = alpha)

	return(output)

}



################################################################################################################
#################### evaluate local polynomial smoother a single x value
################################################################################################################

loc.pol.x <- function(x,data,h,Ker,r)
{

  X <- data$X
  Y <- data$Y
  X.x <- matrix(NA,length(X),0)
  for(j in 0:r)
  {
    X.x <- cbind(X.x,(X-x)^j)
  }
  W.x <- diag( Ker((X-x)/h) )
  a.hat.x <- solve(t(X.x) %*% W.x %*% X.x) %*% t(X.x) %*% W.x %*% Y
  loc.pol.x <- a.hat.x[1]
  return(loc.pol.x)

}

################################################################################################################
#################### get standard errors for local polynomial smoother at single x value
################################################################################################################

loc.pol.x.se <- function(x,data,h,Ker,r,A=NULL)
{	#A is like the projection matrix for a presmoothing estimator
	if(length(A)==0){A <- diag(length(data$X))}

  X <- data$X
  Y <- data$Y
  X.x <- matrix(NA,length(X),0)
  for(j in 0:r)
  {
    X.x <- cbind(X.x,(X-x)^j)
  }
  W.x <- diag( Ker((X-x)/h) )
  cov.a.hat <- solve(t(X.x) %*% W.x %*% X.x) %*% t(X.x) %*% W.x %*% A %*% t(A) %*% W.x %*% X.x %*% solve(t(X.x) %*% W.x %*% X.x)
  loc.pol.se.x <- sqrt(cov.a.hat[1,1])
  return(loc.pol.se.x)
}


################################################################################################################
#################### get local polynomial fit, standard errors, and t.star for confidence bands
################################################################################################################


loc.pol <- function(x, data, h, Ker, r, A=NULL, boot.band = FALSE, e.for.boot.band=NULL, alpha = .05, B = 1000)
{
	if(length(A)==0){A <- diag(length(data$X))}

	# browser()

	n <- length(data$X)

	f.loc.pol <- sapply(data$X,FUN=loc.pol.x,data=data,h=h,Ker=Ker,r=r)
	f.loc.pol.x <- sapply(x,FUN=loc.pol.x,data=data,h=h,Ker=Ker,r=r)
	f.loc.pol.x.se <- sapply(x,FUN=loc.pol.x.se,data=data,h=h,Ker=Ker,r=r,A=A)

	t.star <- NULL

	if( boot.band == TRUE)
	{
		if(length(e.for.boot.band) == 0)
		{

			e.for.boot.band <- data$Y - f.loc.pol

		}

		wboot.e <- matrix( rnorm(n*B,0,sqrt(rep(e.for.boot.band^2,B))),nrow=n)

		delta.wboot.x <- matrix(NA,length(x),B)

		for(b in 1:B)
		{

			delta.wboot.x[,b] <- sapply(x, FUN = loc.pol.x, data = list( X = data$X, Y = as.numeric(A %*% wboot.e[,b]) ), h = h , Ker = Ker, r = r ) / f.loc.pol.x.se

		}

		delta.wboot.sup.over.x <- apply(abs(delta.wboot.x),2,max)

		t.star <- quantile(abs(delta.wboot.sup.over.x),1-alpha)

	}

	output <- list(	f.loc.pol = f.loc.pol,
					f.loc.pol.x = f.loc.pol.x,
					f.loc.pol.x.se = f.loc.pol.x.se,
					x = x,
					t.star = t.star,
					alpha = alpha)

}

################################################################################################################
#################### helper function for crossvalidation for Meier et al estimator
################################################################################################################

Meier_et_al_2009_fold_MSEP <- function(Xtrain,Ytrain,Xtest,Ytest,df,lambda.1.seqs,lambda.2.seq)
{

	n.lambda.1 <- nrow(lambda.1.seqs)
	n.lambda.2 <- length(lambda.2.seq)

	n <- dim(Xtrain)[1]
	q <- dim(Xtrain)[2]

	H <- array(0,dim = c(n,df,q))
	Htest <- matrix(NA,nrow(Xtest),df*q)
	W <- matrix(0,df,df)
	H.tilde <- H.large <- array(0,dim=c(n,q*df,n.lambda.2))
	R.backtransform <- array(0,dim=c(q*df,q*df,n.lambda.2))

	lambda.2.seq <- seq(0,3,length=n.lambda.2)

	all.knots.incl.boundary <- matrix(NA,df-2+6,q)

	# build the design matrices to feed to the grplasso function and transform them to impose the smoothness penalty
	for(j in 1:q)
	{

		knots <- quantile(Xtrain[,j],seq(0,1,length=df-2))
		boundary.knots <- range(knots)
		knots.incl.boundary <- sort(c(rep(boundary.knots,3),knots))

		all.knots.incl.boundary[,j] <- knots.incl.boundary

		H[,,j] <- spline.des(Xtrain[,j], knots = knots.incl.boundary,outer.ok=TRUE,derivs=rep(0,length(Xtrain[,j])))$design
		dsq_bspline.mat <- spline.des(knots, knots = knots.incl.boundary,outer.ok=TRUE,derivs=rep(2,length(knots)))$design

		Htest[,((j-1)*df+1):(df*j)] <- spline.des(Xtest[,j], knots = knots.incl.boundary,outer.ok=TRUE,derivs=rep(0,length(Xtest[,j])))$design

		for(k in 1:df)
			for(l in 1:df)
				{

					pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
					h <- diff(knots)
					W[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(knots)])*h)  # sum of trapezoidal areas.

				}

		for( l in 1:n.lambda.2 )
		{

			M <- 1/n * t(H[,,j]) %*% H[,,j] + lambda.2.seq[l]^2*W
			R <- chol(M)
			H.tilde[,((j-1)*df+1):(df*j),l] <- H[,,j] %*% solve(R)
			R.backtransform[((j-1)*df+1):(df*j),((j-1)*df+1):(df*j),l] <- solve(R)

		}

	}

	grps <- c(NA,c(1:q) %x% rep(1,df))
	Ytest.hat <- array(NA,dim=c(nrow(Xtest),n.lambda.1,n.lambda.2))
	MSEP <- matrix(NA,n.lambda.1,n.lambda.2)

	for( l in 1:n.lambda.2 )
	{

		model <- grplasso(	y = Ytrain,
							x = cbind(rep(1,n),H.tilde[,,l]),
							index = grps,
							lambda = lambda.1.seqs[1:n.lambda.1,l],
							model = LinReg(),
		                  	center = FALSE,
		                    standardize = FALSE,
		                    control = grpl.control(update.hess = "lambda", trace = 0))

		beta.tilde <- model$coef[-1,]

		Ytest.hat <- mean(Ytrain) + Htest %*% R.backtransform[,,l] %*% beta.tilde
		MSEP[,l] <- apply( (Ytest.hat - Ytest %*% matrix(1,1,n.lambda.1))^2,2,mean)

	}

	return(MSEP)

}

################################################################################################################
#################### get cv choices of lambda.1 and lambda.2 for Meier et al estimator
################################################################################################################

Meier_et_al_2009_cv <- function(X,Y,df,n.lambda.1,n.lambda.2,K,x,fois,plot=FALSE,fois.true.x=NA)
{

	n <- dim(X)[1]
	q <- dim(X)[2]

	H <- array(0,dim = c(n,df,q))
	W <- matrix(0,df,df)
	H.tilde <- H.large <- array(0,dim=c(n,q*df,n.lambda.2))
	R.backtransform <- array(0,dim=c(q*df,q*df,n.lambda.2))

	lambda.2.seq <- seq(0,2,length=n.lambda.2)

	all.knots.incl.boundary <- matrix(NA,df-2+6,q)

	# build the design matrices to feed to the grplasso function and transform them to impose the smoothness penalty
	for(j in 1:q)
	{

		knots <- quantile(X[,j],seq(0,1,length=df-2))
		boundary.knots <- range(knots)
		knots.incl.boundary <- sort(c(rep(boundary.knots,3),knots))

		all.knots.incl.boundary[,j] <- knots.incl.boundary

		H[,,j] <- spline.des(X[,j], knots = knots.incl.boundary,outer.ok=TRUE,derivs=rep(0,length(X[,j])))$design
		dsq_bspline.mat <- spline.des(knots, knots = knots.incl.boundary,outer.ok=TRUE,derivs=rep(2,length(knots)))$design

		for(k in 1:df)
			for(l in 1:df)
				{

					pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
					h <- diff(knots)
					W[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(knots)])*h)  # sum of trapezoidal areas.

				}

		for( l in 1:n.lambda.2 )
		{

			M <- 1/n * t(H[,,j]) %*% H[,,j] + lambda.2.seq[l]^2*W
			R <- chol(M)
			H.tilde[,((j-1)*df+1):(df*j),l] <- H[,,j] %*% solve(R)
			R.backtransform[((j-1)*df+1):(df*j),((j-1)*df+1):(df*j),l] <- solve(R)

		}

	}

	grps <- c(NA,c(1:q) %x% rep(1,df))
	beta.tilde <- array(NA,dim=c(df*q,n.lambda.1,n.lambda.2))
	lambda.1.seqs <- matrix(NA,n.lambda.1,n.lambda.2)

	for( l in 1:n.lambda.2 )
	{

		lambda.max <- lambdamax(y = Y,
								x = cbind(rep(1,n),H.tilde[,,l]),
								index = grps,
								model = LinReg(),
			                  	center = FALSE,
			                    standardize = FALSE)

		lambda.exp <- exp(seq(log(lambda.max),0,length = floor(1.1*n.lambda.1))) - 1
		lambda.1.seq <- lambda.exp * lambda.max/max(lambda.exp)

		model <- grplasso(	y = Y,
							x = cbind(rep(1,n),H.tilde[,,l]),
							index = grps,
							lambda = lambda.1.seq[1:n.lambda.1],
							model = LinReg(),
		                  	center = FALSE,
		                    standardize = FALSE,
		                    control = grpl.control(update.hess = "lambda", trace = 0))

		beta.tilde[,,l] <- model$coef[-1,]

		lambda.1.seqs[,l] <- lambda.1.seq[1:n.lambda.1]

	}

	# Estimator fitted on full data set at many values of tuning parameters. Now do crossvalidation.

	indscramble <- sample(1:n,n)
	h <- 1
	n.lambda.1 <- nrow(lambda.1.seqs)
	n.lambda.2 <- length(lambda.2.seq)

	fold.MSEP <- array(NA,dim=c(n.lambda.1,n.lambda.2,K))

	for(k in 1:K)
	{

		fold.ind <- indscramble[h:(h + (floor(n/K)+1)*(k<=(n%%K)) + floor(n/K)*(k>(n%%K)) - 1 )]
		h <- h + (floor(n/K)+1)*(k<=(n%%K)) + floor(n/K)*(k>(n%%K))
		n.fold <- length(fold.ind)

		Xtrain <- X[-fold.ind,]
		Xtest <- X[fold.ind,]
		Ytrain <- Y[-fold.ind]
		Ytest <- Y[fold.ind]


		fold.MSEP[,,k]<- Meier_et_al_2009_fold_MSEP(Xtrain,Ytrain,Xtest,Ytest,df,(K-1)/K*lambda.1.seqs,lambda.2.seq)

		print(k)

	}

	mean.MSEP <- apply(fold.MSEP,c(1,2),mean)
	min.cv <- which( mean.MSEP == min(mean.MSEP), arr.ind = TRUE )

	lambda.1.cv.ind <- min.cv[1]
	lambda.2.cv.ind <- min.cv[2]

	cv.lambda.1 <- lambda.1.seqs[lambda.1.cv.ind,lambda.2.cv.ind]
	cv.lambda.2 <- lambda.2.seq[lambda.2.cv.ind]

	R.backtransform.cv <- R.backtransform[,,lambda.2.cv.ind]
	beta.tilde.cv <- beta.tilde[,lambda.1.cv.ind,lambda.2.cv.ind]

	fois.Meier.x <- matrix(NA,length(x),length(fois))

	for(j in 1:length(fois))
	{
		bspline.mat.x <- spline.des(x = x, knots = all.knots.incl.boundary[,j],outer.ok=TRUE,derivs=rep(0,length(x)))$design
		fois.Meier.x[,j] <- mean(Y) + (bspline.mat.x %*% R.backtransform.cv[((j-1)*df+1):(df*j),((j-1)*df+1):(df*j)]) %*% beta.tilde.cv[((j-1)*df+1):(df*j)]
	}



	if(plot==TRUE)
	{

		par(mfrow=c(2,ceiling(length(fois)/2)),mar=c(5.1, 4.1, 1.1, 1.1))

		for(k in 1:length(fois))
		{

			j <- fois[k]

			plot(NA,ylim=c(-5,5),xlim=c(-2.5,2.5),ylab=paste("f_",j,sep=""),xlab=paste("x_",j,sep=""))

			if(is.na(fois.true.x)==FALSE) { lines(fois.true.x[,j]~x) }

			bspline.mat.x <- spline.des(x = x, knots = all.knots.incl.boundary[,j],outer.ok=TRUE,derivs=rep(0,length(x)))$design

			for(l in 1:n.lambda.2)
				for(m in 1:n.lambda.1)
				{

					f.j.hat.x <- mean(Y) + (bspline.mat.x %*% R.backtransform[((j-1)*df+1):(df*j),((j-1)*df+1):(df*j),l]) %*% beta.tilde[((j-1)*df+1):(df*j),m,l]

					lines(f.j.hat.x~x,col=rgb(1,0,0,(1:m)/m))

				}

			lines(fois.Meier.x[,j]~x,col=rgb(0,0,1,1),lwd=1.5)

		}

	}

	output <- list(	cv.lambda.1 = cv.lambda.1,
					cv.lambda.2 = cv.lambda.2,
					fois.Meier.x = fois.Meier.x,
					x = x)

}

################################################################################################################
#################### compute Meier et al estimator at a single value of lambda.1 and lambda.2
################################################################################################################

Meier_et_al_2009 <- function(X,Y,df,lambda.1,lambda.2,x,fois)
{

	n <- dim(X)[1]
	q <- dim(X)[2]

	H <- array(0,dim = c(n,df,q))
	W <- matrix(0,df,df)
	H.tilde <- H.large <- matrix(0,n,q*df)
	R.backtransform <- matrix(0,q*df,q*df)

	all.knots.incl.boundary <- matrix(NA,df-2+6,q)

	# build the design matrices to feed to the grplasso function and transform them to impose the smoothness penalty
	for(j in 1:q)
	{

		knots <- quantile(X[,j],seq(0,1,length=df-2))
		boundary.knots <- range(knots)
		knots.incl.boundary <- sort(c(rep(boundary.knots,3),knots))

		all.knots.incl.boundary[,j] <- knots.incl.boundary

		H[,,j] <- spline.des(X[,j], knots = knots.incl.boundary,outer.ok=TRUE,derivs=rep(0,length(X[,j])))$design
		dsq_bspline.mat <- spline.des(knots, knots = knots.incl.boundary,outer.ok=TRUE,derivs=rep(2,length(knots)))$design

		for(k in 1:df)
			for(l in 1:df)
				{

					pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
					h <- diff(knots)
					W[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(knots)])*h)  # sum of trapezoidal areas.

				}


			M <- 1/n * t(H[,,j]) %*% H[,,j] + lambda.2^2*W
			R <- chol(M)
			H.tilde[,((j-1)*df+1):(df*j)] <- H[,,j] %*% solve(R)
			R.backtransform[((j-1)*df+1):(df*j),((j-1)*df+1):(df*j)] <- solve(R)

	}

	grps <- c(NA,c(1:q) %x% rep(1,df))

	model <- grplasso(	y = Y,
						x = cbind(rep(1,n),H.tilde),
						index = grps,
						lambda = lambda.1,
						model = LinReg(),
	                  	center = FALSE,
	                    standardize = FALSE,
	                    control = grpl.control(update.hess = "lambda", trace = 0))

	beta.tilde <- model$coef[-1,]

	fois.Meier.x <- matrix(NA,length(x),length(fois))
	for(k in 1:length(fois))
	{
		j <- fois[k]
		bspline.mat.x <- spline.des(x = x, knots = all.knots.incl.boundary[,j],outer.ok=TRUE,derivs=rep(0,length(x)))$design
		fois.Meier.x[,k] <- mean(Y) + (bspline.mat.x %*% R.backtransform[((j-1)*df+1):(df*j),((j-1)*df+1):(df*j)]) %*% beta.tilde[((j-1)*df+1):(df*j)]
	}

	output <- list( lambda.1 = lambda.1,
					lambda.2 = lambda.2,
					fois.Meier.x = fois.Meier.x,
					x = x
					)

}






