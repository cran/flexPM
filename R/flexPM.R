BS <- function (x, df, knots, degree, Boundary.knots = range(x), deriv = 0){
	ord <- 1 + (degree <- as.integer(degree))
	if(ord <= 1){stop("'degree' must be integer >= 1")}
	if(missing(knots)){
		nIknots <- df - ord + 1
		if(nIknots < 0){
			nIknots <- 0
			stop("'df' cannot be smaller than 'degree'")
		}
		knots <- if(nIknots > 0){
			knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
			2)[-c(1, nIknots + 2)]
			stats::quantile(x, knots)
		}
	}
	Aknots <- sort(c(rep(Boundary.knots, ord), knots))
	basis <- splines::spline.des(Aknots, x, ord, derivs = rep(deriv, length(x)))$design
	basis <- basis[, -1L, drop = FALSE]
	a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots, 
		Boundary.knots = Boundary.knots)
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("BS", "basis", "matrix")
	basis
}

flexPM <- function(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6){

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	fit <- flexPM.internal(mf = mf, cl = cl, df = df, 
		degree = degree, knots = knots, maxit = maxit, tol = tol, type = "u")
	fit
}

cflexPM <- function(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6){

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")

	fit <- flexPM.internal(mf = mf, cl = cl, df = df, 
		degree = degree, maxit = maxit, tol = tol, type = "cens")
	fit
}

ctflexPM <- function(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6){

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	fit <- flexPM.internal(mf = mf, cl = cl, df = df, 
		degree = degree, knots = knots, maxit = maxit, tol = tol, type = "trunc")
	fit
}


flexPM.internal <- function(mf,cl, df, degree, knots, maxit, tol, type){

	mt <- attr(mf, "terms")
	yzd <- cbind(model.response(mf))
	n <- nrow(yzd)
	if(type == "u"){y <- yzd; z <- rep.int(-Inf,n); d <- rep.int(1, n)}
	else if(type == "cens"){
		y <- yzd[,1]
		z <- rep.int(-Inf,n)
		d <- (if(ncol(yzd) == 1) rep.int(1, n) else yzd[,2])
	}
	else{	
		y <- yzd[,1]
		z <- yzd[,2]
		z[z < min(y)] <- -Inf
		d <- (if(ncol(yzd) == 2) rep.int(1, n) else yzd[,3])
	}

	x <- model.matrix(mt,mf)
	weights <- model.weights(mf)
	if(!is.null(weights)){weights <- weights/sum(weights)*length(y)}

	model <- L.estimator(y,d,x,weights, df, degree, knots, 
		z = z, maxit = maxit, tol = tol, type = type)
	fit <- list(converged = model$converged, n.it = model$n.it, 
		n = model$n, n.free.par = model$n.free.par, 
		loglik = model$loglik, AIC = model$AIC, BIC = model$BIC
	)
	fit$mf <- mf
	attr(fit$mf, "m") <- model$m
	attr(fit$mf, "s") <- model$s
	attr(fit$mf, "u") <- model$u
	attr(fit$mf, "xystats") <- model$xystats
	fit$call <- cl
	class(fit) <- "flexPM"
	if(!fit$converged){warning("optimization algorithm did not converge")}

	fit
}


L.estimator <- function(y,d,x,weights, df,degree,knots, theta, z, maxit = 100, tol = 1e-6, type){

	if(type == "trunc"){
		if(any(y < z)){stop("y cannot be < z")}
		loglik <- (if(is.null(weights)) L.loglik.cens.trunc else L.wloglik.cens.trunc)
	}
	else if(type == "cens"){loglik <- (if(is.null(weights)) L.loglik.cens else L.wloglik.cens)}
	else{loglik <- (if(is.null(weights)) L.loglik else L.wloglik)}

	# x

	v <- qr((x <- cbind(1,x)))
	sel <- v$pivot[1:v$rank]
	x <- x[,sel, drop = FALSE]
	if((q <- ncol(x)) > 2){
		pcx <- stats::prcomp(x[,-1], center = TRUE, scale = TRUE)
		x[,-1] <- pcx$x
	}
	SdX <- apply(x,2,sd); Mx <- colMeans(x)
	Mx[1] <- 0; SdX[1] <- 1
	x <- scale(x, center = Mx, scale = SdX)
	

	# y, z

	min.y <- min(y)
	max.y <- max(y)
	notrunc <- which(z <= min.y)
	z[notrunc] <- min.y

	y <- (y - min.y)/(max.y - min.y)*10
	z <- (z - min.y)/(max.y - min.y)*10
	if(!missing(knots) && !is.null(knots)){
		knots <- (sort(knots) - min.y)/(max.y - min.y)*10
		if(any(knots <= 0 | knots >= 10)){stop("'knots' must be in the range of the response variable")}
	}
	d1 <- which(d == 1)
	dd <- 1 + d

	# u(.)

	By <- BS(y, df = df, knots = knots, degree = degree, Boundary.knots = c(-1e-5, 10 + 1e-5))
	B1y <- BS(y, df = df, degree = degree, Boundary.knots = c(-1e-5, 10 + 1e-5), deriv = 1, knots = attr(By, "knots"))
	Bz <- BS(z, df = df, degree = degree, Boundary.knots = c(-1e-5, 10 + 1e-5), knots = attr(By, "knots"))

	if((k <- ncol(By)) > 1){
		for(j in (k - 1):1){
			By[,j] <- By[,j] + By[,j + 1]
			B1y[,j] <- B1y[,j] + B1y[,j + 1]
			Bz[,j] <- Bz[,j] + Bz[,j + 1]
		}
	}
	B1y <- pmax(B1y,0)

	# Starting points ##############################################################################

	if(missing(theta)){
		m0 <- lm(rowSums(By) ~ -1 + x, weights = weights)
		mu <- m0$coef
		sigma <- sd(m0$residuals)*sqrt(3)/pi
		theta <- c(mu, log(sigma), rep(0, q - 1), rep(0, k - 1))
	}

	### global optimization ######################################################################

	if(missing(maxit)){maxit <- 10*length(theta)}
	opt <- flexPM.newton(c(theta), loglik, maxit = maxit, tol = tol,
		By = By, B1y = B1y, x = x, d = d, d1 = d1, dd = dd, weights = weights, 
		Bz = Bz, notrunc = notrunc)

	phi <- opt$estimate
	theta <- (if(k > 1) c(exp(phi[(2*q + 1):length(phi)]), 1) else 1)
	mu <- phi[1:q]; sigma <- phi[(q + 1):(2*q)]

	# xy stats

	xystats <- list(
		min.y = min.y, max.y = max.y, mu = mu, sigma = sigma, 
		sel = sel, rot = (if(q > 2) pcx$rotation else NULL), 
		center0 = (if(q > 2) pcx$center else NULL), scale0 = (if(q > 2) pcx$scale else NULL),
		center = attr(x, "scaled:center"), scale = attr(x, "scaled:scale")
	)

	### output

	ll <- -opt$minimum
	n <- length(y)
	r <- length(phi)
	uy <- c(By%*%cbind(theta))
	u <- stats::splinefun(y,uy, method = "monoH.FC")	

	out <- list(
		converged = (opt$n.it < maxit), n.it = opt$n.it, 
		m = c(x%*%cbind(mu)), s = c(exp(x%*%cbind(sigma))) + 0.001,
		u = u, loglik = ll, n = n, n.free.par = r,
		AIC = 2*(r - ll), BIC = -2*ll + r*log(n),
		xystats = xystats
	)
	out
}

# predict CDF, PDF, and survival (type = "CDF"), p-th quantile(s) (type = "QF"), 
# or simulate data with p ~ U(0,1) (type = "sim").
# 'tol' is used to compute quantiles.


predict.flexPM <- function(object, type = c("CDF", "QF", "sim"), newdata, p, ...){
	
	if(is.na(match(type <- type[1], c("CDF", "QF", "sim")))){stop("invalid 'type'")}
	
	mf <- object$mf
	mu <- attr(mf, "m")
	sigma <- attr(mf, "s")
	u <- attr(mf, "u")
	s <- attr(mf, "xystats")
	
	mt <- terms(mf)
	miss <- attr(mf, "na.action")
	nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	xlev <- .getXlevels(mt, mf)

	predCDF <- function(y,u, mu, sigma){
		z1 <- (u(y) - mu)/sigma
		z2 <- log(1 + exp(-z1))
		if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
		-z2
	}
	predQ <- function(u, ya, yb, Fa, Fb, p, mu, sigma){

		noa <- (Fa > p)
		nob <- (Fb < p)
		delta <- 2.5
		while(any(noa) | any(nob)){
			ya[noa] <- ya[noa] - delta
			yb[nob] <- yb[nob] + delta
			Fa[noa] <- predCDF(ya[noa], u, mu[noa], sigma[noa])
			Fb[nob] <- predCDF(yb[nob], u, mu[nob], sigma[nob])
			noa <- (Fa > p)
			nob <- (Fb < p)
			ya[nob] <- yb[nob]
			yb[noa] <- ya[noa]
			delta <- delta*2
		}
		n.it <- ceiling((log(max(yb - ya)) - log(1e-2))/log(2))

		for(i in 1:n.it){
			yc <- (ya + yb)*0.5
			Fc <- predCDF(yc, u, mu, sigma)
			w <- (Fc < p)
			ya[w] <- yc[w]
			yb[!w] <- yc[!w]
		}
		yc <- (ya + yb)*0.5
	}

	if(!missing(newdata)){

		if(type == "CDF"){
			yn <- (if(object$call[[1]] == "flexPM") colnames(mf)[1] else colnames(mf[[1]])[1])
			if(is.na(match(yn, colnames(newdata))))
				{stop("for 'type = CDF', 'newdata' must contain the y-variable")}
			cY <- colnames(Y <- cbind(model.response(mf)))
			if(ncol(Y) > 1){for(j in cY[-1]){newdata[,j] <- 1}}
		}
		else{mt <- delete.response(mt)}
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
			{stop("'newdata' must contain all x-variables")}


		mf <- model.frame(mt, data = newdata, xlev = xlev)
		if(nrow(mf) == 0){
			nr <- nrow(newdata)
			if(type == "CDF"){
				out <- data.frame(matrix(NA,nr,3))
				colnames(out) <- c("log.f", "log.F", "log.S")
				rownames(out) <- rownames(newdata)
			}
			else if(type == "QF"){
				out <- data.frame(matrix(NA,nr,length(p)))
				colnames(out) <- paste("p",p, sep = "")
				rownames(out) <- rownames(newdata)
			}
			else{out <- rep.int(NA, nr)}
			return(out)
		}
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])

		x <- model.matrix(mt, mf)
		x <- cbind(1,x)
		x <- x[,s$sel, drop = FALSE]
		if((q <- ncol(x)) > 2)
			{x[,-1] <- scale(x[,-1, drop = FALSE], scale = s$scale0, center = s$center0)%*%s$rot}
		x <- scale(x, center = s$center, scale = s$scale)
		mu <- c(x%*%cbind(s$mu))
		sigma <- c(exp(x%*%cbind(s$sigma))) + 0.001
	}

	if(type == "CDF"){
		y <- cbind(model.response(mf))[,1]
		y <- (y - s$min.y)/(s$max.y - s$min.y)*10
		z1 <- (u(y) - mu)/sigma
		z2 <- log(1 + exp(-z1))
		z3 <- log(1 + exp(z1))
		if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
		if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
		log.F <- log.S <- log.f <- NULL
		log.F[nomiss] <- -z2
		log.S[nomiss] <- -z3
		log.f[nomiss] <- -z1 - log(sigma) - 2*z2 + log(u(y, deriv = 1)) + log(10) - log(s$max.y - s$min.y)
		log.F[miss] <- log.S[miss] <- log.f[miss] <- NA
		out <- data.frame(log.f = log.f, log.F = log.F, log.S = log.S)
		rownames(out)[nomiss] <- rownames(mf)
		if(!is.null(miss)){rownames(out)[miss] <- names(miss)}
		return(out)
	}
	else{
		ya.0 <- yb.0 <- rep.int(5, length(mu))
		Fa.0 <- Fb.0 <- predCDF(ya.0, u, mu, sigma)

		if(type == "QF"){
			if(missing(p)){stop("please indicate 'p'")}
			if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
			P <- log(p)
			out <- NULL
			for(j in 1:length(p)){
				Qp <- NULL
				Qp[nomiss] <- predQ(u, ya.0, yb.0, Fa.0, Fb.0, P[j], mu, sigma)
				Qp[miss] <- NA
				out <- cbind(out, Qp*(s$max.y - s$min.y)/10 + s$min.y)
			}
			colnames(out) <- paste("p", p, sep = "")
			out <- as.data.frame(out)
			rownames(out)[nomiss] <- rownames(mf)
			if(!is.null(miss)){rownames(out)[miss] <- names(miss)}
			return(out)
		}
		else{
			p <- log(runif(length(mu)))
			Qp <- NULL
			Qp[nomiss] <- predQ(u, ya.0, yb.0, Fa.0, Fb.0, p, mu, sigma)
			Qp[miss] <- NA
			Qp <- Qp*(s$max.y - s$min.y)/10 + s$min.y
			return(Qp)
		}
	}
}

# print method

print.flexPM <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of obs: ", paste(deparse(round(x$n)), sep = " ", collapse = " "), "\n", sep = "")
  cat("n. of free parameters: ", paste(deparse(round(x$n.free.par)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Converged: ", paste(deparse(x$converged), sep = " ", collapse = " "), "\n", sep = "")
  cat("N. of iterations: ", paste(deparse(round(x$n.it)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Log-likelihood: ", paste(deparse(round(x$loglik,1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("AIC: ", paste(deparse(round(x$AIC,1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("BIC: ", paste(deparse(round(x$BIC,1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Use predict() to obtain predictions from the fitted model")
  cat("\n")
}

# Some u(y) is assumed to be logistic
# This is the log-likelihood function, weigthed/unweighted and truncated/non truncated.
# If deriv = 0, only returns the log-likelihood
# If deriv = 1, returns the N*q matrix of first derivatives
# If deriv = 2, returns the gradient and the hessian


L.wloglik <- function(theta, By, B1y, x,d,d1,dd,weights, Bz, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:(ncol(By) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(ncol(By) > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))

	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	loglik <- sum(ll*weights)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*2 - 1
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a*weights/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - 1)*weights)) # 
	if(ncol(By) > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(1/u1y) - By.h*(a/sigma) #
		d.theta <- ColSums(d.theta*weights)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*2*weights/sigma))
	H_sigma <- t(x)%*%(x*(weights*((v^2)*z1*(2*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - 1)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*2 - a/sigma)*weights))

	if(ncol(By) > 1){
		H_theta <- (diag(d.theta, ncol(By) - 1, ncol(By) - 1) 
			- t(B1y.h)%*%(B1y.h*(weights/u1y^2)) - t(By.h)%*%(By.h*(weights*2/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*2*weights/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(weights*v*(a/sigma + z1*2*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}

L.loglik <- function(theta, By, B1y, x,d,d1,dd,weights, Bz, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:(ncol(By) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(ncol(By) > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))


	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	loglik <- sum(ll)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*2 - 1
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - 1)))
	if(ncol(By) > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma)
		d.theta <- ColSums(d.theta)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*2/sigma))
	H_sigma <- t(x)%*%(x*(((v^2)*z1*(2*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - 1)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*2 - a/sigma)))

	if(ncol(By) > 1){
		H_theta <- (diag(d.theta, ncol(By) - 1, ncol(By) - 1) 
			- t(B1y.h)%*%(B1y.h*(1/u1y^2)) - t(By.h)%*%(By.h*(2/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*2/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(v*(a/sigma + z1*2*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}


L.wloglik.cens <- function(theta, By, B1y, x,d,d1,dd,weights, Bz, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:(ncol(By) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(ncol(By) > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))

	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll*weights)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*dd - d
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a*weights/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d)*weights))
	if(ncol(By) > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma)
		d.theta <- ColSums(d.theta*weights)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*dd*weights/sigma))
	H_sigma <- t(x)%*%(x*(weights*((v^2)*z1*(dd*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - d)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - a/sigma)*weights))

	if(ncol(By) > 1){
		H_theta <- (diag(d.theta, ncol(By) - 1, ncol(By) - 1) 
			- t(B1y.h)%*%(B1y.h*(weights*d/u1y^2)) - t(By.h)%*%(By.h*(weights*dd/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*dd*weights/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(weights*v*(a/sigma + z1*dd*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}

L.loglik.cens <- function(theta, By, B1y, x,d,d1,dd,weights, Bz, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:(ncol(By) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(ncol(By) > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))


	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*dd - d
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d)))
	if(ncol(By) > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma)
		d.theta <- ColSums(d.theta)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*dd/sigma))
	H_sigma <- t(x)%*%(x*(((v^2)*z1*(dd*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - d)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - a/sigma)))

	if(ncol(By) > 1){
		H_theta <- (diag(d.theta, ncol(By) - 1, ncol(By) - 1) 
			- t(B1y.h)%*%(B1y.h*(d/u1y^2)) - t(By.h)%*%(By.h*(dd/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*dd/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(v*(a/sigma + z1*dd*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}


L.wloglik.cens.trunc <- function(theta, By, B1y, x,d,d1,dd,weights, Bz, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)

	h <- 1:(ncol(By) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(ncol(By) > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))
	uz <- c(Bz%*%cbind(theta))

	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))	
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}

	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)

	# truncation
	
	zz <- (uz - mu)/sigma
	log.Fz <- -log(1 + exp(-zz))
	log.Sz <- -log(1 + exp(zz))
	if(any(outFz <- (log.Fz == -Inf))){log.Fz[outFz] <- zz[outFz]}
	if(any(outSz <- (log.Sz == -Inf))){log.Sz[outSz] <- -zz[outSz]}
	log.Sz[notrunc] <- 0
	log.Fz[notrunc] <- -Inf

		
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll*weights) - sum(log.Sz*weights)

	###
	
	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*dd - d
	az <- exp(log.Fz)
	a_az <- a - az
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*((a_az)*weights/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d - az*zz)*weights))
	if(ncol(By) > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		Bz.h <- Bz[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma) + Bz.h*(az/sigma)
		d.theta <- ColSums(d.theta*weights)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
	dFz_dg <- exp(log.Sz + log.Fz - log.sigma)
	dFz_dmu <- -dFz_dg
	dFz_dsigma <- dFz_dmu*zz


	H_mu <- t(x)%*%(x*((dF_dmu*dd - dFz_dmu)*weights/sigma))
	H_sigma <- t(x)%*%(x*(weights*((v^2)*(z1*(dd*sigma*dF_dsigma - a) - zz*(sigma*dFz_dsigma - az)) + v*(0.001/sigma*(a*z1 - d - az*zz)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - dFz_dsigma - a_az/sigma)*weights))

	if(ncol(By) > 1){
		H_theta <- (diag(d.theta, ncol(By) - 1, ncol(By) - 1) 
			- t(B1y.h)%*%(B1y.h*(weights*d/u1y^2)) 
			- t(By.h)%*%(By.h*(weights*dd/sigma*dF_dg)) 
			+ t(Bz.h)%*%(Bz.h*(weights/sigma*dFz_dg)))
		H_mu_theta <- t(x)%*%((By.h*(dF_dg*dd) - Bz.h*dFz_dg)*weights/sigma)
		H_sigma_theta <- t(x)%*%(By.h*(weights*v*(a/sigma + z1*dd*dF_dg)) - Bz.h*(weights*v*(az/sigma + zz*dFz_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}


L.loglik.cens.trunc <- function(theta, By, B1y, x,d,d1,dd,weights, Bz, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)

	h <- 1:(ncol(By) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(ncol(By) > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))
	uz <- c(Bz%*%cbind(theta))

	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))	
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}

	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)

	# truncation
	
	zz <- (uz - mu)/sigma
	log.Fz <- -log(1 + exp(-zz))
	log.Sz <- -log(1 + exp(zz))
	if(any(outFz <- (log.Fz == -Inf))){log.Fz[outFz] <- zz[outFz]}
	if(any(outSz <- (log.Sz == -Inf))){log.Sz[outSz] <- -zz[outSz]}
	log.Sz[notrunc] <- 0
	log.Fz[notrunc] <- -Inf

		
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll) - sum(log.Sz)

	###
	
	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*dd - d
	az <- exp(log.Fz)
	a_az <- a - az
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*((a_az)/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d - az*zz)))
	if(ncol(By) > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		Bz.h <- Bz[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma) + Bz.h*(az/sigma)
		d.theta <- ColSums(d.theta)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
	dFz_dg <- exp(log.Sz + log.Fz - log.sigma)
	dFz_dmu <- -dFz_dg
	dFz_dsigma <- dFz_dmu*zz


	H_mu <- t(x)%*%(x*((dF_dmu*dd - dFz_dmu)/sigma))
	H_sigma <- t(x)%*%(x*(((v^2)*(z1*(dd*sigma*dF_dsigma - a) - zz*(sigma*dFz_dsigma - az)) + v*(0.001/sigma*(a*z1 - d - az*zz)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - dFz_dsigma - a_az/sigma)))

	if(ncol(By) > 1){
		H_theta <- (diag(d.theta, ncol(By) - 1, ncol(By) - 1) 
			- t(B1y.h)%*%(B1y.h*(d/u1y^2)) 
			- t(By.h)%*%(By.h*(dd/sigma*dF_dg)) 
			+ t(Bz.h)%*%(Bz.h*(1/sigma*dFz_dg)))
		H_mu_theta <- t(x)%*%((By.h*(dF_dg*dd) - Bz.h*dFz_dg)/sigma)
		H_sigma_theta <- t(x)%*%(By.h*(v*(a/sigma + z1*dd*dF_dg)) - Bz.h*(v*(az/sigma + zz*dFz_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}


flexPM.newton <- function(start, f, tol = 1e-5, maxit = 200, ...){

	f0 <- f(start, ..., deriv = 2)
	g <- attr(f0, "gradient")
	h <- attr(f0, "hessian")
	conv <- FALSE
	eps <- 1
	alg <- "nr"

	for(i in 1:maxit){

		if(conv | max(abs(g)) < tol){break}

		####
		
		H1 <- try(chol(h), silent = TRUE)
		if(class(H1) != "try-error"){
			if(alg == "gs"){alg <- "nr"; eps <- 1}
			delta <- chol2inv(H1)%*%g
		}
		else{
			if(alg == "nr"){alg <- "gs"; eps <- 1}
			delta <- g
		}

		####

		f1 <- Inf
		while(f1 > f0){
			new.start <- start - delta*eps
			if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
			f1 <- f(new.start, ..., deriv = 0)
			eps <- eps*0.5
			if(is.na(f1)){f1 <- Inf}
		}

		if(conv | f0 - f1 < tol){break}
		f1 <- f(new.start, ..., deriv = 2)
		g1 <- attr(f1, "gradient")
		h1 <- attr(f1, "hessian")
		start <- new.start; f0 <- f1; g <- g1; h <- h1
		eps <- min(eps*10,1)
	}

	list(estimate = start, n.it = i, minimum = as.numeric(f0), gradient = g, hessian = h)
}




























