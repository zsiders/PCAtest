#' PCoAtest: Statistical Significance of PCoA
#'
#'   PCoAtest uses random permutations to build null distributions for several statistics of a PCA analysis: Psi (Vieira 2012), Phi (Gleason and Staelin 1975), the rank-of-roots (ter Braak 1988). Comparing these distributions with the observed values of the statistics, the function tests: (1) the hypothesis that there is more correlational structure among the observed variables than expected by random chance, (2) the statistical significance of each PC, and (3) the contribution of each observed variable to each significant PC. The function also calculates the sampling variance around mean observed statistics based on bootstrap replicates.
#'
#' @param x A matrix or dataframe with variables in the columns and the observations in the rows.
#' @param dfun Distance measurement function for PCoA
#' @param nperm Number of random permutations to build null distributions of the statistics.
#' @param nboot Number of bootstrap replicates to build 95%-confidence intervals of the observed statistics.
#' @param alpha Nominal alpha level for statistical tests.
#' @param plot A logical specifying whether to plot the null distributions, observed statistics, and 95%-confidence intervals of statistics based on random permutation and bootstrap resampling.
#' @param corr correction for PCoA, "none", "lingoes", "cailliez"
#' @param print A logical (TRUE/FALSE) indicating whether to print benchmarks to console
#'
#' @details
#' PCoAtest uses the function ape::pcoa to run a PCoA. PCoAtest plots three types of graphs in a single page: (1) a histogram showing the null distribution and the observed value of the Psi statistic, (2) a histogram showing the null distribution and the observed value of the Phi statistic, and (3) a bar plot of the percentage of explained variance of each PC1, PC2, ..., etc., showing the sampling variance based on bootstrap replicates and random permutations with 95%-confidence intervals. If the PCA is not significant, based on the Psi and Phi testing results, only histograms (1) and (2) are shown.
#'
#' @return
#' An object of class “list” with the following elements:
#' \describe{
#'   \item{emp_psi}{The observed Psi statistic.}
#'
#'   \item{emp_phi}{The observed Phi statistic.}
#'
#'   \item{null_psi}{The null distribution of Psi values.}
#'
#'   \item{null_phi}{The null distribution of Phi values.}
#' 
#' 	 \item{eigenprob}{Probabilities of empirical eigenvalues}
#' 
#' 	 \item{sigpcs}{Logical vector of significant PCs}
#'
#'   \item{percvar_emp}{The percentage of variance explained by each PC based on the observed data.}
#'
#'   \item{percvar_boot}{The percentage of variance explained by each PC based on the bootstrapped data.}
#'
#'   \item{percvar_perm}{The percentage of data explained by each PC based on the permuted data.}
#'
#'}
#'
#' @author
#' Arley Camargo (original author)
#' Zachary A. Siders (PCoA & parallel, zsiders@ufl.edu)
#'
#' @references
#' \itemize{
#'   \item Gleason, T. C. and Staelin R. (1975) A proposal for handling missing data. Psychometrika, 40, 229–252.
#'   \item Ringnér, M. (2008) What is principal component analysis? Nature Biotechnology, 26, 303–304.
#'   \item ter Braak, C. F. J. (1990) Update notes: CANOCO (version 3.1). Agricultural Mattematic Group, Report LWA-88-02, Wagningen, Netherlands.
#'   \item Vieira, V. M. N. C. S. (2012) Permutation tests to estimate significances on Principal Components Analysis. Computational Ecology and Software, 2, 103–123.
#'   \item Wong, M. K. L. and Carmona, C. P. (2021) Including intraspecific trait variability to avoid distortion of functional diversity and ecological inference: Lessons from natural assemblages. Methods in Ecology and Evolution. \url{https://doi.org/10.1111/2041- 210X.13568}.
#'}
#' @export
#'
#' @examples
#'#PCA analysis of five uncorrelated (r=0) variables
#'library(MASS)
#'mu <- rep(0,5)
#'Sigma <- matrix(c(rep(c(1,0,0,0,0,0),4),1),5)
#'ex0 <- mvrnorm(100, mu = mu, Sigma = Sigma )
#'result<-PCoAtest(ex0, dfun = dist, 100, 100, 0.05, plot=TRUE)
#'
#'#PCA analysis of five correlated (r=0.25) variables
#'Sigma <- matrix(c(rep(c(1,0.25,0.25,0.25,0.25,0.25),4),1),5)
#'ex025 <- mvrnorm(100, mu = mu, Sigma = Sigma )
#'result<-PCoAtest(ex025, dfun = dist, 100, 100, 0.05, plot=TRUE)
#'
#'#PCA analysis of five correlated (r=0.5) variables
#'Sigma <- matrix(c(rep(c(1,0.5,0.5,0.5,0.5,0.5),4),1),5)
#'ex05 <- mvrnorm(100, mu = mu, Sigma = Sigma )
#'result<-PCoAtest(ex05, dfun = dist, 100, 100, 0.05, plot=TRUE)
#'
#'#PCoA analysis of seven morphological variables from 29 ant species (from
#'#Wong and Carmona 2021, https://doi.org/10.1111/2041-210X.13568)
#'data("ants")
#'result<-PCoAtest(ants, dfun = dist, 100, 100, 0.05, plot=TRUE)

PCoAtest <- function(x, dfun = dist, nperm=1000, nboot=1000, alpha=0.05, plot=TRUE, cores.free=4, corr = 'lingoes', print = TRUE) {

	# check dependencies
	requireNamespace("ape", quietly = TRUE)
	requireNamespace("stats", quietly = TRUE)
	requireNamespace("grDevices", quietly = TRUE)
	requireNamespace("graphics", quietly = TRUE)
	requireNamespace("utils", quietly = TRUE)

	#run function
		eig_calc <- function(x = x, dfun = dfun, corr = corr){
			d <- dfun(x)
			pcoa <- ape::pcoa(d, correction=corr)
			eigenvalues <- pcoa$values$Eigenvalues[pcoa$values$Eigenvalues>0]

			if (dim(x)[1] < dim(x)[2]) {
				eigenvalues <- eigenvalues[-length(eigenvalues)]
			}
			pervarobs <- eigenvalues/sum(eigenvalues) * 100

			return(list(eigenval = eigenvalues,
			            pervar = pervarobs))
		}

		psiphi_calc <- function(x, eigenvalues){
			Psi <- 0
			for (i in 1:length(eigenvalues)) {
			  Psi <-  Psi + (eigenvalues[i]-1)^2
			}

			Phi <- 0
			for (i in 1:length(eigenvalues)) {
			  Phi <-  Phi + eigenvalues[i]^2
			}

			Phi <- sqrt((Phi - dim(x)[2]) / (dim(x)[2] * (dim(x)[2]-1)))

			return(list(Psi = Psi, Phi = Phi))
		}

	# empirical eigenvalues, loadings, Psi, and Phi
		
		eigobs <- eig_calc(x, dist, corr)
		eigenobs <- eigobs$eigenval
		pervarobs <- eigobs$pervar
		ppobs <- psiphi_calc(x, eigenobs)
		Psiobs <- ppobs$Psi
		Phiobs <- ppobs$Phi

	# bootstraped data for confidence intervals of empirical eigenvalues, index loadings,
	# and correlations

		if(print) cat("\nSampling bootstrap replicates... Please wait\n")

		library(doSNOW)
		cl <- makeCluster(parallel::detectCores()-cores.free)
		registerDoSNOW(cl)
		pb <- txtProgressBar(max = nboot, style = 3)
		progress <- function(n) setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
		pervarboot <- foreach(i = 1:nboot,.combine = rbind, .options.snow = opts) %dopar%{
		    bootdata <- x[sample(nrow(x),size=dim(x)[1],replace=TRUE),]
				eigboot <- eig_calc(bootdata, dfun, corr)

				return(eigboot$pervar)
		}
		close(pb)
		stopCluster(cl)


		# confidence intervals of percentage of variation
		confint <- apply(pervarboot,MARGIN=2,FUN=stats::quantile, probs=c(0.025,0.975))


	# null distributions based on randomizations of Psi, Phi, eigenvalues, percentage of variation,
	# and index loadings


		if(print)cat("\nSampling random permutations... Please wait\n")
		library(doSNOW)
		cl <- makeCluster(parallel::detectCores()-cores.free)
		registerDoSNOW(cl)
		pb <- txtProgressBar(max = nperm, style = 3)
		progress <- function(n) setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
		permres <- foreach(i = 1:nperm, .options.snow = opts) %dopar%{
		    repvalue <- 0
				perm <- apply(x, MARGIN=2, FUN=sample)
				eigperm <- eig_calc(perm, dfun, corr)
				ppperm <- psiphi_calc(perm, eigperm$eigenval)
				return(list(eigenval = eigperm$eigenval,
				            pervar = eigperm$pervar,
				            psi = ppperm$Psi,
				            phi = ppperm$Phi))
		}
		close(pb)
		stopCluster(cl)

		eigenrand <- do.call('rbind', lapply(permres, function(x) x$eigenval))
		pervarperm <- do.call('rbind', lapply(permres, function(x) x$pervar))
		Psiperm <- do.call('c', lapply(permres, function(x) x$psi))
		Phiperm <- do.call('c', lapply(permres, function(x) x$phi))


		confintperm <- apply (pervarperm, MARGIN=2, FUN=stats::quantile, probs=c(0.025,0.975))

		# comparing empirical Psi, Phi and eigenvalues with their null distributions to calculate
		# p-values

		Psiprob <- length (which (Psiperm>Psiobs)) / nperm
		Phiprob <- length (which (Phiperm>Phiobs)) / nperm

		eigenprob <- rep(NA,length(eigenobs))
		for (k in 1:length( eigenobs)) {
			eigenprob[k] <- length (which (eigenrand[,k] > eigenobs[k])) / nperm
		}

		if(Psiprob < alpha & Phiprob < alpha) { # test PC axes if both Psi and Phi are significant

		# find out which PCs are significant

			sigaxes <- 0
			for (i in 1:length(eigenprob)) {
		  	if (eigenprob[i] < alpha) {
		  		sigaxes <- sigaxes + 1
		  	}
			}
		}

	# screen output

		if(print){
			cat(paste("\n", "========================================================", sep=""),
			paste("Test of PCA significance: ", dim(x)[2], " variables, ", dim(x)[1], " observations", sep=""),
			paste(nboot, " bootstrap replicates, ", nperm, " random permutations", sep=""),
			paste("========================================================", sep=""), sep="\n")


			cat(paste("\n", "Empirical Psi = ", format(round(Psiobs,4),nsmall=4), ", Max null Psi = ", format(round(max(Psiperm),4),nsmall=4), ", Min null Psi = ", format(round(min(Psiperm),4),nsmall=4), ", p-value = ", Psiprob, sep=""),
				paste("Empirical Phi = ", format(round(Phiobs,4),nsmall=4),", Max null Phi = ", format(round(max(Phiperm),4),nsmall=4),", Min null Phi = ", format(round(min(Phiperm),4),nsmall=4),", p-value = ", Phiprob, sep=""), sep="\n")

			if (Psiprob >= alpha & Phiprob >= alpha) { # if both Psi and Phi are not significant
			  cat(paste ("\n", "PCA is not significant!", sep=""))
			}

			if (Psiprob < alpha & Phiprob < alpha) { # test PC axes if both Psi and Phi are significant

				for (i in 1:length(eigenobs)) {
					cat(paste ("\n", "Empirical eigenvalue #", i, " = ", round(eigenobs[i],5), ", Max null eigenvalue = ", round(max(eigenrand[,i]), 5), ", p-value = ", eigenprob[i], sep=""))
				}

				cat("\n")

				for (i in 1:sigaxes) {
					cat(paste ("\n", "PC ", i, " is significant and accounts for ", round(pervarobs[i], digits=1), "% (95%-CI:",round(confint[1,i],digits=1),"-",round(confint[2,i],digits=1),")"," of the total variation", sep=""))
				}

				if (sigaxes > 1) {
					cat("\n")
					cat(paste ("\n", "The first ", sigaxes, " PC axes are significant and account for ", round(sum(pervarobs[1:sigaxes]), digits=1), "% of the total variation", sep=""))
				}

				cat("\n\n")
			}
		}

	# plots

	if (plot==T) {
		graphics::par(mfrow=c(2,2), mar=c(5, 4, 1, 2) + 0.1)

		# plot of empirical and randomized Psi

		h <- graphics::hist(Psiperm,plot=FALSE)
		h$density = h$counts/sum(h$counts)*100
		plot(h,freq=FALSE, col="gray45", xlab="Psi", xlim=c(0, max(max(Psiperm), Psiobs)), ylab="Percentage of permutations", main="")
		graphics::arrows(x0=Psiobs, y0=max(h$density)/10, y1=0, col="red", lwd=2, length=0.1)
		graphics::legend("topright", legend=c("Null distribution","Empirical value"), fill=c("gray45","red"), bty="n", cex=0.8)

		# plot of empirical and randomized Phi

		h <- graphics::hist(Phiperm,plot=FALSE)
		h$density <- h$counts/sum(h$counts)*100
		plot(h,freq=FALSE, col="gray45", xlab="Phi", xlim=c(0, max(max(Phiperm), Phiobs)), ylab="Percentage of permutations", main="")
		graphics::arrows (x0=Phiobs, y0=max(h$density)/10, y1=0, col="red", lwd=2, length=0.1)

		# plot of bootstrapped and randomized percentage of variation for all PC's

		if (Psiprob < alpha & Phiprob < alpha) { # test PC axes if both Psi and Phi are significant

			plot(pervarobs, ylab="Percentage of total variation", xlab="PC", bty="n", ylim=c(0, max(confint)), type="b", pch=19, lty="dashed", col="red", xaxt = "n")
		  graphics::axis(1, at = 1:length(eigenobs))
			graphics::lines (apply(pervarperm,MARGIN=2,FUN=mean), type="b", pch=19, lty="dashed", col="gray45")
			suppressWarnings(graphics::arrows (x0=c(1:length(eigenobs)), y0=confint[2,], y1=confint[1,], code=3, angle=90, length=0.05, col="red"))
			suppressWarnings(graphics::arrows (x0=c(1:length(eigenobs)), y0=confintperm[2,], y1=confintperm[1,], code=3, angle=90, length=0.05, col="gray45"))
		}
	}

	results <- list()

	results[["emp_psi"]] <- Psiobs
	results[["emp_phi"]] <- Phiobs
	results[["null_psi"]] <- Psiperm
	results[["null_phi"]] <- Phiperm
	results[["eigenprobs"]] <- eigenprob
	results[["sigpcs"]] <- eigenprob < alpha

	if (Psiprob < alpha & Phiprob < alpha) { # test PC axes if both Psi and Phi are significant

		results[["percvar_emp"]] <- pervarobs
		results[["percvar_boot"]] <- pervarboot
		results[["percvar_perm"]] <- pervarperm

	}

	return (results)

}
