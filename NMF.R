<source-file filename="NMF.R" display="link" linktext="NMF and data-processing routines in R">

# NMF Library routines

pnor <- function(x) {
	x/sum(x)
}


PlotM<- function(Z,heatc=12) {
m<-ncol(Z);
n<-nrow(Z);
MM<-matrix(0,m,n);
MM[1:m,]<-t(Z[,1:m][n:1,]);
image(MM,axes=FALSE,col=heat.colors(heatc));
}

# Normalize a matrix
#
# Args:
# 	m: the matix
# 	L: by rows (L=1) or columns (L=2)
# Returns:
# 	The normalized matrix
NormalizeWH <- function(m, L = 2) {
	m2 <- apply(m, L, innp <- function(x) {
			return (x / sqrt(sum(x^2)))
			})
	mtmn <- t(m2) %*% m2
	return(mtmn)
}

# Non-negative matrix factorization
#
# Args:
#		difcon: indicates how difference between the similarity matrice in two iteration.
#too small, in many times it has not convergenced in 2000 steps 
#too large, because it convergenced in 200 steps
#may try to see results on 1e-5,-6,-7,-8,-9,-10, to see if it is sensentive on convergence criteria
NMF <- function(Z, k, method = "lee", eps = 1e-16, stopconv = 40, num.iter = 2000, difcon = 1e-10, ifseed = TRUE, seed = 10) {
	n <- nrow(Z)
	m <- ncol(Z)
	# Connectivity matrix
	cons <- matrix(rep(0, m*m),ncol = m, nrow = m)
	consold <- cons
	inc <- 0

	if (min(Z) < 0) {
		stop("No negative elements allowed")
	}

	# Test if there is void data: all elements in a row are zeros.
	if(nrow(as.matrix(which(rowSums(Z) > 0))) != n) {
		stop("Zero rows are not allowed")
	}

# initiation randomly select matrix values from uniform distribution of [min(Z),max(Z)]
#set.seed()
	if (ifseed == TRUE) {
		set.seed(seed)
	}
	H <- runif(k*m, min = min(Z), max = max(Z))
	H <- matrix(H, k, m)

	if (ifseed == TRUE) {
		set.seed(seed)
	}
	W <- runif(k*n, min = min(Z), max = max(Z))
	W <- matrix(W, n, k)

	for (i in 1:num.iter){
		# Adjust small values to avoid undeflow (i.e., avoid dividing zero)
		if (i %% 10 == 0){
			maxmat <- function(x) {
				max(x, eps)
			}
			H <- apply(H, c(1,2), maxmat)
			W <- apply(W, c(1,2), maxmat)
		}

		WH1 <- W %*% H

		# Update H
		if (method == "brunet") {
			sumW <- matrix(rep(t(colSums(W)), m), k, m)
			H <- H * (t(W) %*% (Z / WH1)) / sumW
                        H[is.na(H)] <- eps
		}
		if (method == "lee") {
			H <- H * ((t(W) %*% Z) / (t(W) %*% W %*% H))
                        H[is.na(H)] <- eps
		}

		WH2 <- W %*% H

		# Update W  
		if (method=="brunet") {
			sumH <- matrix(rep(t(rowSums(H)), n), n, k, byrow = TRUE)
			W <- W * ((Z / WH2) %*% t(H)) / sumH
                        W[is.na(H)] <- eps
		}
		if (method == "lee") {
			W <- W * ((Z %*% t(H)) / (W %*% H %*% t(H)))
                        W[is.na(H)] <- eps
		}

		# Construct connectivity matrix
		cons <- NormalizeWH(H)
		difcontest <- sum((cons - consold)^2) / (m^2)
		consold <- cons
		if (is.nan(difcontest)) {
			break
		}
		if (difcontest < difcon) {
			# Connectivity matrix has not changed
			inc <- inc + 1
		} else {
			inc <- 0
		}
		if (inc > stopconv){
			break
		}
	}

	difdistance <- sum((Z - W %*% H)^2)
	#print(paste("iteration =", i, "difdistance =", difdistance))

	if (i == num.iter) {
		print("NMF has not converged after maximum number of allowed iterations")
	}
	return(list(W = W, H = H))
}

# Find the best optimum solution for NMF give an rank K
#
# Args:
#		Z: the input matrix
#		K: the rank
#		nloop: the iterative steps
# provide one method to select W and H
# WBest and HBest means we select the factorization with the smallest KL
BestKLDiver <- function(Z, K, nloop, nmf.method = "brunet", r.nmf = FALSE, difcon = 1e-10, ifseed = FALSE) {
	n <- nrow(Z)
	m <- ncol(Z)

	Wmatrix0 <- array(0,dim=c(nloop,n,K))
	Hmatrix0 <- array(0,dim=c(nloop,K,m))

	spw0 <- rep(0,n)
	sph0 <- rep(0,m)

	SD <- rep(0,nloop)
	for (j in 1:nloop) {
	print(paste("Computing iteration", j, "of", nloop, "for rank", K))

		if (r.nmf == TRUE){
                        library(NMF) # pre-installed R package "NMF" 
			L <- nmf(Z, K, method = nmf.method);
			#print(summary(L)[5])
			W <- basis(L)
			H <- coef(L)
		} else {
			L <- NMF(Z,K,method=nmf.method,difcon=difcon,ifseed=ifseed,seed=j)
			W <- L$W
			H <- L$H
		}

		SD[j] <- sum((Z - W %*% H)^2)
		Wmatrix0[j,,] <- W
		Hmatrix0[j,,] <- H
	}

	a <- which(SD == min(SD))
	#print(a)
	Wmatrix <- Wmatrix0[a,,]
	Hmatrix <- Hmatrix0[a,,]
	spw <- apply(Wmatrix, 1, sparseness)
	sph <- apply(Hmatrix, 2, sparseness)

	return(list(SD = SD,
		sph = sph,
		spw = spw,
		Wmatrix = Wmatrix,
		Hmatrix = Hmatrix))

}

KLDive <- function(X, WH2){
   	WH2[WH2 == 0] <- 1e-300
		X[which(X == 0)] <- 1e-300
		D <- sum(X * log(X / WH2) - X + WH2)
		return(D)
}


sparseness <- function(x) {
	n <- length(x)
	sparsedegree <- (sqrt(n) - sum(abs(x)) / sqrt(sum(x^2))) / (sqrt(n) - 1)
	return(sparsedegree)
}

# Select taxa similar to a basis pfam
write.taxa.list <- function(X = An, H = H, P = P, labelnames, taxaid, pfams, pfam.file, spw1, ixv = ixv, k, simiv = 0.8, profile, V.m) {
	ix <- which(X[, k] > simiv)
	IX <- sort(X[ix, k], decreasing = TRUE, index.return = TRUE)
	ix2 <- IX$ix
	ix <- ix[ix2]
	M <- P[ix,ixv]
	M <- data.frame(taxaid = taxaid[ix], taxa = pfams[ix], similarity = X[ix, k])
	rownames(M) <- NULL
	write.table(M,con<-pfam.file)

	return(ix = ix)
}

AffineMatrix <- function(HTHN, r = 0.2){
	A <- exp(-((1 - HTHN)^2) / (2 * (r^2)))
	return(A)
}

SpectralReord <- function(HTH, method = "Lap", evrank = 1) {
	n <- nrow(HTH)
	a <- colSums(HTH)
	D <- diag(a)
	DT <- solve(D)^(1/2)
 
	if (method == "Lap") {
		L <- D - HTH
	} else {
		L <- DT %*% HTH %*% DT
	}
	eg <- eigen(L)
	Ix <- sort(eg$vectors[, n - evrank], index.return = TRUE)
	HTH2 <- HTH[Ix$ix, Ix$ix]
	diag(HTH2) <- 1
	return(list(ix = Ix$ix,
		egva = eg$values,
		ev = eg$vectors[, n - evrank],
		egv = eg$vectors,
		HTH2 = HTH2))
}

indexdiff <- function(n) {
	x <- y <- c(1:n)
	B <- matrix(0,n,n)
	for (i in 1:n) {
		for (j in 1:n) {
			B[i, j] <- abs(x[i] - y[j])
		}
	}
	return(B)
}


# Perform a spectral reordering
#
# Args:
# 	H: the matrix to reorder
# 	L: If the matrix is a W matrix, L = 1; if the matrix is an H matrix, L = 2.
spectreorder <- function(H, HTH, L, beta = c(0.01, 0.3), lim = 0.01) {
	
        bseq <- seq(beta[1], beta[2], lim)
	K <- length(bseq)
       
        if (H){
	HTH <- NormalizeWH(H, L)
        }
        else {HTH=HTH}

	diag(HTH) <- 1
	n <- nrow(HTH)

	C <- rep(0, K)
	D <- rep(0, K)
	v <- matrix(0, n, K)
	evmatrix <- matrix(0, n, K)
	for (i in 1:K) {
		HTHN <- AffineMatrix(HTH, r = bseq[i])
		Slap <- SpectralReord(HTHN)
		ixv <- Slap$ix
		ev <- Slap$ev
		HTHr <- HTH[ixv, ixv]
		B <- indexdiff(n)
		C[i] <- sum(HTHr * (B^2))
		D[i] <- sum(HTHN * (B^2))
		v[,i] <- ixv
		evmatrix[, i] <- ev 
	}

	Min <- which.min(C)
	V <- v[, Min]
	ev <- evmatrix[, Min]

	HTHN <- AffineMatrix(HTH, r = bseq[Min])

	return(list(C = C,
		V = V,
		HTH = HTH[V, V],
		D = D,
		HTHN = HTHN,
		v = v,
		ev = ev))
}
rnor <- function(x) {
	x/sqrt(sum(x^2))
}

ConsensusFuzzyH <- function(Z, kstart, kend, nloop, method = "square", Rnmf = TRUE, nmf.method = "brunet", difcon = 1e-10, ifseed = FALSE) {

	n <- nrow(Z);
	m <- ncol(Z);
	# Test for negative values;
	if (min(Z) < 0) {
		stop('Some elements are negative!')
	}

	if (nrow(as.matrix(which(rowSums(Z) > 0))) != n) {
		stop('A row is all zero!')
	}
		
	KL <- rep(0, kend - 1)
	EUD <- rep(0, kend - 1)
	averdiff <- rep(0, kend - 1)

	for (j in kstart:kend) { 
		s <- 0
		t <- 0
		V <- nloop * (nloop - 1) / 2
		difference <- rep(0, V)
		connvec <- array(0, c(nloop, m, m))
		E <- rep(0, nloop)
		D <- rep(0, nloop)
		for (i in 1:nloop) {
			s <- s + 1
	         print(paste("factor =", j, "iteration =", i))

			if (Rnmf == TRUE) {
				if (ifseed == TRUE) {
					set.seed(i)
				}
                                library(NMF)
				L <- nmf(Z, j, method = nmf.method)
				W <- basis(L)
				H <- coef(L)
			} else {
				L <- NMF(Z, j, method = nmf.method,
                                         difcon = difcon, ifseed = ifseed, seed = i)
				W <- L$W
				H <- L$H
			}

			D[i] <- KLDive(Z, W %*% H)	
			E[i] <- sum((Z - W %*% H)^2)
			connh <- NormalizeWH(H)                   
			connvec[i,,] <- connh
                       
			if (s > 1) {
				for (k in 1:(s - 1)){
					if (method == "abs") {
						a <- sum(abs(connh-connvec[k,,]))
					} else if (method == "square") {
						a <- sum((connh - connvec[k, , ])^2)
					}
					t <- t + 1
                                       # print(a)
					difference[t] <- a
				}
			}
		}
               #print(difference)
		KL[j - 1] <- sum(D[i]) / nloop
		EUD[j - 1] <- sum(E[i]) / nloop
		averdiff[j - 1] <- 1 - sum(difference) / (V * (m^2))
		#print(paste("KL =", KL[j - 1], "EUD =", EUD[j - 1], "averdiff =", averdiff[j - 1]))
	}

	return(list(averdiff = averdiff,
				KL = KL,
				EUD = EUD))
}

# Replace the image function
source("http://lalashan.mcmaster.ca/Rstuff/newimage.R")

# a function to plot a matrix by labeling columns and rows.
biplots <- function(OrdZ,pathways,samples,colorset=12,op= par(mar=c(10,2,4,8)), shown=100, LabelAngel=90, colorf=heat.colors,cexaxis=0.7,cexayis=0.7){
#nrowZ: the number of rows of Z for visulization
	nrowZ<-nrow(OrdZ)
        ncolZ <-ncol(OrdZ)
	min <- min(OrdZ)
	max <- max(OrdZ)
       
#rainbow, red background, green blue, not good
#heat.colors, red and yellow,not good
#terrain.colors,  green background, yellow
#topo.colors,  blue background, green
#cm.colors.  light blue background, cym 

	ColorRamp <-colorf(colorset);
	ColorLevels <- seq(min, max, length=length(ColorRamp))

# Reverse Y axis
		reverse <- c(nrowZ:1)
		creverse <- c(ncolZ:1)
		OrdZ <- OrdZ[reverse,creverse]
                        if (!is.logical(pathways)){
	                   yLabels <- pathways
                           yLabels <- yLabels[reverse]}#supposed to be pathways

                      	if (!is.logical(samples)){xLabels <- samples;}

		op<- op
                
                if (!is.logical(samples)){        
		image(1:length(xLabels), 
				1:length(yLabels), 
				t(OrdZ), 
				col=ColorRamp,
				xlab="", ylab="", 
				axes=FALSE, 
				zlim=c(min,max),
				font.axis=2, useRaster=TRUE)}
                if (is.logical(samples)){
                image(          t(OrdZ), 
				col=ColorRamp,
				xlab="", ylab="", 
				axes=FALSE, 
				zlim=c(min,max),
				font.axis=2, useRaster=TRUE)}

                if (is.logical(samples)){} else{
		if (ncolZ<shown){ axis(BELOW<-1, 
				at=1:length(xLabels),las=2, 
				labels= xLabels[creverse],
				cex.axis=cexaxis);
                                } else{}}

                if (is.logical(pathways)){} else{
		if (nrowZ<shown){ 
			axis(LEFT <-4,
					at=1:length(yLabels),
					labels=yLabels, las= HORIZONTAL<-1,
					pos<-1, cex.axis=cexayis);

                                 } else{}}
   	par(op)
		layout(1)

}



##############
#selecting pfams

GO_class <- function(GOC="P",PfamX,x=x)
{       
			#For a given pfam (selected) x, provide a GO annotation if it exists.
			
			#GOC indicates the GO classification, "P" indicates biological process, "F" molecular function, and "C" cellular components. 
         #Pfam2X: Pfams and their GO annotation

			#pfamgolist: pfams ID in GO annotations
         pfamgolist <- as.vector(PfamX[,1]) # pfam ID in GO annotation list 
			
			GO <- as.character("NA")
			listgo <- which((pfamgolist==as.character(x))
			                 &(PfamX[,5]==GOC))
			GOlist <- PfamX[listgo,3]
			GO <- as.character(GOlist[1])
         # print(GOlist)
			if (length(GOlist)>1){
		   	for (i in 2:length(GOlist)){
		   	GO <- paste(as.character(GO),as.character(GOlist[i]),sep=";")
		   	}
			}
			return(GO)
}

rnor <- function(x) {
	x/sqrt(sum(x^2))
}

GO_Assign <- function(pfamv, PfamX)
{
         # assigne GO classification to pfamv ( the selected pfam list)
         if (length(pfamv) > 1){
           GOlistM <- lapply(pfamv,
			 function(x){
		                c(GO_class(GOC="P",PfamX,x=x), 
				 GO_class(GOC="F",PfamX,x=x), 
				 GO_class(GOC="C",PfamX,x=x))
	 		 }
		       )
           GOlistM <- t(matrix(unlist(GOlistM),nrow=3))
           
          }
          else{
           GOlistM <- matrix(c(GO_class(GOC="P",PfamX,x=x), 
		        GO_class(GOC="F",PfamX,x=x), 
		        GO_class(GOC="C",PfamX,x=x)),nrow=1)}

         colnames(GOlistM) <- c("Biological_Process", "Molecular_Function", "Cellular_Components")
         return(GOlistM)
}


Selection_list_p  <- function(X,closedegree=closedegree,nMat=nMat, GO.annotation=TRUE, delNA=TRUE,PfamX=PfamX,
                              pfam_file=pfam_file, pfams=pfams) 
{        
         #select pfam list for all components and write to a .txt file
         ix <- which(X>closedegree) #Select the top specific pfams 
         if (length(ix) >0 ){
			IX <- sort(X[ix],decreasing=TRUE,index.return=TRUE) # reorder them by decreasing order
         ix <- ix[IX$ix]       
         M <- nMat[ix,] # reordered profile

         #print("-------------")
         #print(paste("High similar pfams to "," basis Pfam profile "))
         
			pfamv <- as.vector(pfams[ix,1]) # ID of selected pfams
                        
         if (GO.annotation==TRUE){
		   	GOlistM <- GO_Assign(pfamv, PfamX)
                      if (delNA==TRUE){
		   	#delet "NA" rows
		      ls.v <- apply(GOlistM, 1, function(x){
		      	return(ifelse((is.na(x[1])&is.na(x[2])&is.na(x[3])),0,1))
		      })
                      ls.na <- which(ls.v==1) # at least one GO classification is not "NA"
                      } 

                      if (delNA==FALSE) {ls.na <- 1:nrow(GOlistM)}

		   	M <- data.frame(PFAM_ID= (pfams[ix,1])[ls.na], PFAM_Description=(pfams[ix,3])[ls.na],
                            Closeness_to_BPP=round(X[ix],digits=3)[ls.na], GO =GOlistM[ls.na,]
                           )
			}
		   if (GO.annotation==FALSE){
         	M <- data.frame(PFAM_ID= (pfams[ix,1]), PFAM_Description=(pfams[ix,3]), Closeness_to_BPP=round(X[ix],digits=3)
                           )
			}
		   rownames(M) <- NULL
			#write to pfam_file
			cat("\n\n\n", file=pfam_file)
         cat("\n----------------------",file=pfam_file)
         cat(paste("High similar pfams to"," basis Pfam profile (BPP)"),file=pfam_file)
         cat("-----------------------\n\n\n",file=pfam_file)
         write.table(M,con<-pfam_file)
			Closeness <- round(X[ix],digits=3)
         }

			if (length(ix)==0){
			ix <- 0
			Closeness <- "None items matched"
			}
         return(list(ix = ix,Closeness=Closeness))
}


Selection_list_k  <- function(X, closedegree=closedegree, closeness=closeness, nMat=nMat, delNA=TRUE, csvfile=csvfile, GO.annotation=TRUE, PfamX=PfamX, 
                              pfam_file=pfam_file, pfams=pfams,Horder=Horder,k=k,ntable=ntable) 
{        
         if ( is.matrix(X) ) { Xk <- X[,k]}
         if ( !is.matrix(X) )  {Xk <- X}
         #select pfam list for the component k and write to a .txt file
         ix <- which(Xk>closedegree) #For a given k to select the top specific pfams 
         if (length(ix) >0 ){
         IX <- sort(Xk[ix],decreasing=TRUE,index.return=TRUE) # reorder them by decreasing order
         ix <- ix[IX$ix]       
         M <- nMat[ix,] # reordered profile

         #print("-------------")
         #print(paste("High similar pfams to "," basis Pfam profile ",k))
         
			pfamv <- as.vector(pfams[ix,1]) # ID of selected pfams
                        
         if (GO.annotation==TRUE){
		   	GOlistM <- GO_Assign(pfamv, PfamX)
		   	#delet "NA" rows
                      if (delNA==TRUE){
		      ls.v <- apply(GOlistM, 1, function(x){
		      	return(ifelse((is.na(x[1])&is.na(x[2])&is.na(x[3])),0,1))
		      })
                      ls.na <- which(ls.v==1) # at least one GO classification is not "NA" 
                      } 

                      if (delNA==FALSE) {ls.na <- 1:nrow(GOlistM)}
		      
		   	M <- data.frame(PFAM_ID= (pfams[ix,1])[ls.na], PFAM_Description=(pfams[ix,3])[ls.na],
                            Closeness_to_BPP=round(Xk[ix],digits=3)[ls.na], GO =GOlistM[ls.na,]
                           )
			   if (nrow(M) >ntable){
				M <- M[1:ntable,]
				}
			}
		   if (GO.annotation==FALSE){
         	M <- data.frame(PFAM_ID= (pfams[ix,1]), PFAM_Description=(pfams[ix,3]), Closeness_to_BPP=round(Xk[ix],digits=3)
                           )
			}

		rownames(M) <- NULL

	#write to pfam_file
         if (!is.logical(csvfile)){
            if ( !is.matrix(X) ){
            write.table(M[1:30,], file=paste(csvfile,".tsv",sep=""), sep="\t");
            write.table(M, file=paste(csvfile,".csv",sep=""), sep="\t")}
            if ( is.matrix(X) ){
            write.table(M[1:30,], file=paste(csvfile,k,".tsv",sep=""), sep="\t");
            write.table(M, file=paste(csvfile,k,".csv",sep=""), sep="\t")}
         }

			cat("\n\n\n", file=pfam_file)
         cat("\n----------------------",file=pfam_file)
         cat(paste("High similar pfams to"," basis Pfam profile (BPP)",Horder[k]),file=pfam_file)
         cat("-----------------------\n\n\n",file=pfam_file)
         write.table(M,con<-pfam_file)
			Closeness <- round(Xk[ix],digits=3)
         }

			if (length(ix)==0){
			ix <- 0
			Closeness <- "None items matched"
			}
         return(list(ix = ix,Closeness=Closeness))
}

         
Selection_list <- function(nMat=nMat,Wmatrix=Wmatrix, Hmatrix=Hmatrix, delNA=TRUE,
                            closeness= "similarity", Horder = c(1,2,3,4,5), csvfile = csvfile, 
									 pfam_file =pfam_file, pfams, closedegree = 0.9, GO.annotation=GO.annotation, Pfam2X=Pfam2X,ntable=50)
{
	      #A general function to select pfam list based on how a pfam specific to a component under a given "similarity measure".
          
		   #nMat <- apply(Mat,2,function(x){x/sum(x)}) # normalization of columns
                   #nMat
			profile <- t(apply(as.matrix(nMat),1,rnor)) # normalization of rows
 
			K <- nrow(Hmatrix) 
			H <- apply(Hmatrix[Horder,],1,rnor) # normalization of rows 

			if (closeness=="similarity"){
          An <- profile %*% H  # similarity matrix of pfams and basis Pfam profiles
			}
         
			if (closeness=="correlation"){ 
          An <- cor(t(profile),H) # correlation matrix of pfams and basis Pfam profiles
			}

			if (closeness=="specificity"){
			#Specificity method based on W only
			 sparseness<-function(x){
             n <- length(x);
             sparsedegree<-(sqrt(n)-sum(abs(x))/sqrt(t(x)%*%x))/(sqrt(n)-1)
             return(sparsedegree)}
          An <- apply(Wmatrix,1,sparseness)
			}

         if (closeness=="entropy"){
		   #Gene score method based on W only 
			# See: Sparse Non-negative Matrix Factorizations via Alternating Non-negativity-constrained Least Squares for Microarray Data Analysis, Hyunsoo Kim and Haesun Park, Bioinformatics, 23-12:1495-1502, 2007. 
			  prob <- Wmatrix/apply(Wmatrix,1,sum)
           An <- 1 + apply(prob*log2(prob),1,sum)/log2(ncol(Wmatrix))	
			}

			if (closeness=="projection"){
			# selecting pfams based on Wmatrix
			An <- Wmatrix
                        #/apply(Wmatrix,1,sum)
			}
			
         # This is the give file names to print the annotated list
         pfam_file <- file(pfam_file,"w") #

         An <- as.matrix(An)

         if (ncol(An)>1){
         Closeness <- NULL
         indk <- NULL
         ind1 <- NULL
         lix0 <- 1
         for (k in 1:K){
           ixS <- Selection_list_k(X=An, closeness=closeness, closedegree=closedegree,nMat=nMat, pfam_file=pfam_file, csvfile=csvfile, PfamX=Pfam2X, delNA=delNA,
			                          pfams=pfams,Horder=Horder,k=k, GO.annotation=GO.annotation,ntable=ntable)
           ix <- ixS$ix
           S <-  ixS$Closeness
			  lix <- length(ix)

           if (lix <= ntable){
            if (lix >0)
            {
             indk[lix0:(lix0+lix-1)] <- k
             ind1[lix0:(lix0+lix-1)] <- ix
             Closeness[lix0:(lix0+lix-1)] <- S
             lix0=lix0+lix;
            }
            else{
             lix=1
             indk[lix0:(lix0+lix-1)] <- k
             ind1[lix0:(lix0+lix-1)] <- ix
             Closeness[lix0:(lix0+lix-1)] <- S
             lix0=lix0+lix;
            }
           }
           
           if (length(ix) >ntable){
			    indk[lix0:(lix0+ntable-1)] <- k
             ind1[lix0:(lix0+ntable-1)] <- ix[1:ntable]
             Closeness[lix0:(lix0+ntable-1)] <- S[1:ntable]
             lix0=lix0+ntable;
			  }

        }

        close(pfam_file) # close the file
        ind <- cbind(indk,ind1,Closeness)}
        
        if (ncol(An)==1) {
        
		  ixP <- Selection_list_p(X=An,closedegree=closedegree,nMat=nMat, delNA=delNA,
		                   GO.annotation=TRUE, PfamX=Pfam2X, pfam_file=pfam_file, pfams=pfams) 

		    
                    lc <- which(An > closedegree)
                    ixlc <- sort(An[lc],decreasing=TRUE,index.return=TRUE)
	            llc <- length(lc)
		    if ((llc >0)&(llc<=ntable)){
		      ind <- data.frame(1:llc,(ixlc$ix), Closeness=(ixlc$x))
			 }
			 if (llc>ntable){
			   ind <- data.frame(1:(ntable*5),(ixlc$ix)[1:(ntable*5)], Closeness=(ixlc$x)[1:(ntable*5)])
			 }
			 if (length(lc)==0){
		 	 ind <- c(0,0,0)
			 }
		   }
        
		  return(list(ind=ind))
}

# to plot a colorgrid

ColorScale <- function(colmatrix=colmatrix,levelnum=levelnum,colorpanel=c("red","yellow")){
   # generate color level firstly, it has levelnum levels between min(colmatrix) and max(colmatrix)
   mincol<-min(colmatrix)
   maxcol<-max(colmatrix)
   if (maxcol == mincol)
     maxcol<-mincol+1
   col <- colorRampPalette(colorpanel)(levelnum)
   ColorLevels <- seq(mincol, maxcol, length=levelnum)
   return(list(ColorLevels=ColorLevels,col=col))
}


cr = colorRamp(heat.colors(100))

valCol = function(v){
	rgb(cr(v), maxColorValue=255)
}

rangescale <- function(v) {(v-min(v))/(max(v)-min(v))}

areascale <- function(v) {sqrt(v/max(v))/2} # Max radius = 1/2

wrapLabels <- function(labels, width) {
  rejoin <- function(lines) { paste(lines,collapse="\n") }
  return( sapply( strwrap(labels, width=width, simplify=FALSE), rejoin ) )
}


# circlePlot is used by bubblePlot

circlePlot <- function(x, y, sample, bg,levelnum=5,HTH,label,cexaxis=0.3){
 
       layout(matrix(data=c(1,2,3,4), nrow=2,ncol=2,byrow=TRUE),
		 widths=c(2.5,2),heights=c(1,2))
       
       nrowZ<-nrow(HTH)
       min <- min(HTH)
       max <- max(HTH)
 
       yLabels <- label
       xLabels <- label
 
     #  jet.colors <- colorRampPalette(
     #             Col <-  c( "blue","white", "red"))
      
       ColorRamp <-heat.colors(100);
       ColorLevels <- seq(min, max, length=length(ColorRamp))
 
        # Reverse Y axis
        reverse <- c(nrowZ:1)
        yLabels <- yLabels[reverse]
        HTH <- HTH[reverse,]
       
        par(mar=c(0,1.5,13,6.5))
        image(ColorLevels,1,matrix(data=ColorLevels,
		  nrow=length(ColorLevels),ncol=1),col=ColorRamp,xlab="",
		  ylab="",yaxt="n",cex.axis=1.3, useRaster=TRUE)
        
		  # a fake figure
		 # par(mar=c(5,5,5,5))
		  plot(0,0,col="white",
		  xlab="",
		  ylab="",
		  axes=FALSE,
		  plot=FALSE,
		  main="",
		  yaxt="n",
		  xaxt="n")

        par(mar=c(10,1.5,3,6.5))
        image(1:length(xLabels), 1:length(yLabels), t(HTH), col=ColorRamp, xlab="", ylab="", axes=FALSE, zlim=c(min,max),font.axis=2, useRaster=TRUE)
        axis(BELOW<-1, at=1:length(xLabels),las=2, labels=xLabels, cex.axis=cexaxis)
        if (nrowZ<100){ axis(LEFT <-4, at=1:length(yLabels),labels=yLabels, las= HORIZONTAL<-1,pos<-3, cex.axis=cexaxis)} else{}
 
       par(mar = c(9,1,2.5,7), bg = "white")
        a = 0
		  plot(
        c(a, a+max(as.numeric(x))),
        c(a, a+max(as.numeric(y))),
        xlab="", ylab="", asp=1, type='n', axes=FALSE,
        #main="Environment parameter"
    )

       axis(BELOW<-1, at=1:6,las=2, labels=c("Salinity","Sample Depth","Chlorophyll","Temperature","Water Depth","Insolation"), cex.axis=0.4)

        symbols(as.numeric(x), as.numeric(y),
        squares=2*areascale(sample), 
        bg=bg,
        inches=FALSE,
        add=TRUE
    )

   text(11.5,42,"General legend",cex=0.5)    
   legend(8,40,c("Very high","high","Medium","Low","Not avalaible"),
          pch=15,col=c("green","white","yellow","red","black"),cex=0.5,box.col="white")

   text(12.5,22,"Water depth legend",cex=0.5)    
   legend(8,20,c(">4000m","2000-4000m","900-2000m","~150m","~70m","~1m"),
          pch=15,col=c("green","blue","lightskyblue","white","yellow","red"),cex=0.5,box.col="white")
}

#this function is stolen from http://www.r-bloggers.com/ggheat-a-ggplot2-style-heatmap-function/
## m=matrix(data=sample(rnorm(100,mean=0,sd=2)), ncol=10)
## this function makes a graphically appealing heatmap (no dendrogram) using ggplot
## whilst it contains fewer options than gplots::heatmap.2 I prefer its style and flexibility
 
ggheat=function(m, rescaling='none', clustering='none', labCol=T, labRow=T, border=FALSE, 
heatscale= c(low = "red", mid="yellow", high = "white"), midpoint=0.5) 
{
  ## the function can be be viewed as a two step process
  ## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
  ## using simple options or by a user supplied function
  ## 2. with the now resahped data the plot, the chosen labels and plot style are built
 
  require(reshape)
  require(ggplot2)
 
  ## you can either scale by row or column not both! 
  ## if you wish to scale by both or use a differen scale method then simply supply a scale
  ## function instead NB scale is a base funct
 
  if(is.function(rescaling))
  { 
    m=rescaling(m)
  } 
  else 
  {
    if(rescaling=='column') 
      m=scale(m, center=T)
    if(rescaling=='row') 
      m=t(scale(t(m),center=T))
  }
 
  ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
  ## if you want a different distance/cluster method-- or to cluster and then scale
  ## then you can supply a custom function 
 
  if(is.function(clustering)) 
  {
    m=clustering(m)
  }else
  {
  if(clustering=='row')
    m=m[hclust(dist(m))$order, ]
  if(clustering=='column')  
    m=m[,hclust(dist(t(m)))$order]
  if(clustering=='both')
    m=m[hclust(dist(m))$order ,hclust(dist(t(m)))$order]
  }
	## this is just reshaping into a ggplot format matrix and making a ggplot layer
 
  rows=dim(m)[1]
  cols=dim(m)[2]
  melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows) ,melt(m))
	g=ggplot(data=melt.m)
 
  ## add the heat tiles with or without a white border for clarity
 
  if(border==TRUE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='white')
  if(border==FALSE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
 
  ## add axis labels either supplied or from the colnames rownames of the matrix
 
  if(labCol==T) 
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
  if(labCol==F) 
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
 
  if(labRow==T) 
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
	if(labRow==F) 
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
 
  ## get rid of grey panel background and gridlines
 
  g2=g2+opts(panel.grid.minor=theme_line(colour=NA), panel.grid.major=theme_line(colour=NA),
  panel.background=theme_rect(fill=NA, colour=NA))
 
  ## finally add the fill colour ramp of your choice (default is red to yellow to white)-- and return
  return(g2+scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3], midpoint=midpoint))
 
}
 
  ## NB because ggheat returns an ordinary ggplot you can add ggplot tweaks post-production e.g. 
  ## data(mtcars)
  ## x= as.matrix(mtcars)
  ## ggheat(x, clustCol=T)+ opts(panel.background=theme_rect(fill='pink'))

# Replace the image function
source("http://lalashan.mcmaster.ca/Rstuff/newimage.R")
biplotz <- function(OrdZ,l,pathways,samples,colorset=100,par1=par(mar=c(7,5,5,5)), font.lab=2,  xlab="",ylab="", colorbar= TRUE,
                    modu2x= 6,colorf=heat.colors,cexaxis=0.6, axislab=TRUE, par2=par(mar=c(15,1,10,8)),
						  cexaxisx=0.7,cexaxisy=0.4,
                    boundary=FALSE){

      #l is the scale for visualization
      min <- min(OrdZ)
      max <- max(OrdZ)
      OrdZ <- OrdZ^l
      nrowZ<-nrow(OrdZ)
      ncolZ<- ncol(OrdZ)
		wrapLabels <- function(labels, width) {
        rejoin <- function(lines) { paste(lines,collapse="\n") }
        return( sapply( strwrap(labels, width=width, simplify=FALSE), rejoin ) )
      }
      
      yLabels <-pathways # or wrapLabels(pathways,50)
		xLabels <- samples;
 
      ColorRamp <-colorf(colorset);
      ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- c(nrowZ:1)
 creverse <- c(ncolZ:1)
 yLabels <- yLabels[reverse]
 xLabels <- xLabels[creverse]
 OrdZ <- OrdZ[reverse,creverse]

if (colorbar==TRUE){
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(3.5,1), heights=c(1,2.3))
}
op1<- par1

 image(1:length(xLabels), 1:length(yLabels), t(OrdZ), col=ColorRamp, xlab=xlab, ylab=ylab, axes=FALSE, zlim=c(min,max),font.axis=2, font.lab=font.lab, useRaster=TRUE)

 if (axislab==TRUE){
 axis(BELOW<-1, at=1:length(xLabels),las=2, labels=xLabels, cex.axis=cexaxisx)
 if (nrowZ<200){ axis(LEFT <-4, at=1:length(yLabels),labels=yLabels, las= HORIZONTAL<-1,pos<-3, cex.axis=cexaxisy)} else{}}

   if (is.list(boundary)){
      moduid=1:length(boundary)
      lb <- length(boundary)
     for (ls in 1:lb){
       a <- boundary[[ls]]
      if (ls==1){
       text(a[1]-2.2,a[3]+modu2x*3,paste("M_",moduid[ls],sep=""),cex=1.2,font=2)
     }

       if(ls==lb){
       #text(a[1],a[3]+2,paste("M_",ls,sep=""))
        text(a[2]+3.5,a[3]+modu2x*4,paste("M_",moduid[ls],sep=""),cex=1.2,font=2)
       }

       if ((ls>1) &(ls<lb)){
      text(a[1]-2.2,a[3]+modu2x*4,paste("M_",moduid[ls],sep=""),cex=1.2,font=2)

       }
      
      tick = 0.5
      segments(a[1]+tick,a[3]+tick,a[2]+tick,a[3]+tick,lwd=1.5,lty=2)
      segments(a[1]+tick,a[4]+tick,a[2]+tick,a[4]+tick,lwd=1.5,lty=2)
      segments(a[2]+tick,a[4]+tick,a[2]+tick,a[3]+tick,lwd=1.5,lty=2)
      segments(a[1]+tick,a[3]+tick,a[1]+tick,a[4]+tick,lwd=1.5,lty=2)
    }}

      par(op1)
if (colorbar==TRUE){
  op2<-par2
  # Color Scale
 image(1, ColorLevels, matrix(data=(ColorLevels^l),  ncol=length(ColorLevels),nrow=1),col=ColorRamp,font.lab=font.lab, xlab="",ylab="",xaxt="n", useRaster=TRUE)
 layout(1)}

}

##The function for plot multiple ggplot figures in one plot.
## see here for reference: http://gettinggeneticsdone.blogspot.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE, widths=c(1,1), heights=c(1,1)) {
 dots <- list(...)
 n <- length(dots)
 if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
 if(is.null(nrow)) { nrow = ceiling(n/ncol)}
 if(is.null(ncol)) { ncol = ceiling(n/nrow)}
        ## NOTE see n2mfrow in grDevices for possible alternative
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow,ncol,widths=widths, heights=heights) ) )
 ii.p <- 1
 for(ii.row in seq(1, nrow)){
 ii.table.row <- ii.row 
 if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
  for(ii.col in seq(1, ncol)){
   ii.table <- ii.p
   if(ii.p > n) break
   print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
   ii.p <- ii.p + 1
  }
 }
}



#Function for translate Longitude and Latitude to degree
#some weired symbols in GOS metadata are moved 
Dimension2Degree <- function(D){

#R:
#negative West South
#positive North East

 South <- grep("S",D)
 West <- grep("W",D)

 n <- length(D) 
 x1 <- sub("∞","_",D)
 x2 <- sub("\'","_",x1)
 x <- sub("[A-Z]\"","",x2)
 z <- strsplit(x,"_")
 d <- rep(0,n)

 for (i in 1:n){
 zvector <- z[[i]]
 zvector <- as.numeric(zvector)
 if (length(zvector)<2) {zvector[2:3]<-c(0,0)}
 if (length(zvector)<3) {zvector[3]<-0}
 zvector[which(is.na(zvector)==TRUE)]<-0
 d[i] <- zvector[1]+zvector[2]/60 + zvector[3]/3600
 }

 d[South] <- -d[South]
 d[West] <- -d[West]
 d
}



## plot similarity matrix on world map
plotmap <- function(HTH,ixv, meta=Gos_meta_filtered_rast,pcex,ppch,
                    orientation,addline=TRUE,Th){

#required packages
library(maps)
library(mapproj)

#input HTH: similarity matrix; ixv, order; Gos_meta_filtered_rast: metadata from GOS
 
m <- map("world",plot=FALSE)
 
map('world',proj='azequalarea',col=4,orient=orientation)
#map.grid(m,col=2)
 
 
#read Longitude and modify a small issue
x0 <- as.vector(meta$Longitude)
x0[length(x0)] <- paste(as.character(x0[length(x0)]),"\"",sep="")
 
#read Latitude and modify a small issue
y0 <- as.vector(meta$Latitude)
y0[length(x0)] <- paste(as.character(y0[length(y0)]),"\"",sep="")
y0[length(x0)-1] <- paste(as.character(y0[length(y0)-1]),"\"",sep="")
 
# x longitude; y latitude
# x and y are two vectors giving longitude and latitude coordinates of points on the earth
 
x <- Dimension2Degree(x0)
y <- Dimension2Degree(y0)
 
#make a map project, changing degree to coordinates
Pfam_proj <- mapproject(x[ixv],y[ixv],proj='azequalarea',orient=orientation)
 
n <- nrow(HTH)
 
#make a colormap for nodes along the similarity matrix
cr = colorRamp(c("red", "yellow","white"))
valCol = function(v){
    rgb(cr(v), maxColorValue=255)
}
 
color <- valCol(0:(n-1)/n)
 
#color[3:n] <- valCol(Ev[3:n])
#color[1:2] <- "green"
 
X <- Pfam_proj$x
Y <- Pfam_proj$y
 
GeoLocation <- meta$Geographic.Location
GeoLocation <- GeoLocation[ixv]
Geo <- unique(GeoLocation)
 
S <- unlist(lapply(Geo,function(x){
    l<- which(GeoLocation==x)
    l[c(1)]
  }
))
 
GeoLocation[-S] <- NA
 
GeoLocation
 
#GeoLocation <- (meta$Sample.Name)[ixv]
 
points(Pfam_proj,col=color,pch=ppch,cex=pcex)
points(Pfam_proj,col="black",pch=1,cex=pcex)
 
 for (i in 1:n){
 
   if (addline){
    for (j in 1:n){
      S_h <- HTH[i,j]  
      lwd <- S_h
		if (S_h > Th){
      segments(X[i],Y[i],X[j],Y[j],lwd=lwd,col="green")
		}
      }
     }
}
 
#plot sites
points(Pfam_proj,col=color,pch=ppch,cex=pcex)
points(Pfam_proj,col="black",pch=1,cex=pcex)
 
for (i in 1:n){
 text(X[i]+0.3,Y[i],GeoLocation[i],cex=0.7)
}
 
}


## This function is analogous to the function circlePlot in this page, but EnvGrid only plot the environmental parameters.
EnvGrid <- function(x, y, sample, mar=c(5,4,4,2)+0.1, bg,levelnum=5, cexaxis=0.3){
 
        par(mar=mar)
        a = 0
		  plot(
        c(a, a+max(as.numeric(x))),
        c(a, a+max(as.numeric(y))),
        xlab="", ylab="", asp=1, type='n', axes=FALSE,
        #main="Environment parameter"
    )
 
       axis(BELOW<-1, at=1:6,las=2, font=2,labels=c("Salinity","Sample Depth","Chlorophyll","Temperature","Insolation", "Water Depth"),col.axis="grey50", cex.axis=cexaxis)
 
        symbols(as.numeric(x), as.numeric(y),
        squares=2*areascale(sample), 
        bg=bg,
        inches=FALSE,
        add=TRUE
    )
 
   text(11.5,42,"General",cex=0.8)    
  # text(11.5,42,"General",cex=0.8)    
   legend(8,40,c("High","Medium","Low","NA"), text.col="gray50",
          pch=15,col=c("white","yellow","red","black"),cex=0.8,box.col="white")
   legend(8,40,c("High","Medium","Low","NA"), text.col="gray50",
          pch=0,col=c("black","black","black","black"),cex=0.8,box.col="white")

   text(12.5,22,"Water depth",cex=0.8)
 #  text(12.5,22,"Water depth",cex=0.8)     
   legend(8,20,c(">4000m","2000-4000m","900-2000m","100-200m","20-100m","0-20m"),text.col="gray50",
          pch=15,col=c("darkblue","blue3","blue2","blue","lightblue3","lightblue"),cex=0.8,box.col="white")
   legend(8,20,c(">4000m","2000-4000m","900-2000m","100-200m","20-100m","0-20m"),text.col="gray50",
          pch=0,col=c("black","black","black","black","black","black"),cex=0.8,box.col="white")
}

                #lu <- which(grepl("unknown",M$PFAM_Description)==TRUE)
                #if (length(lu)<1){MList <- M} else {
                #unknownList <- M[lu,]
                #knownList <- M[-lu,]

                #rownames(unknownList) <- NULL
                #rownames(knownList) <- NULL
                #nk <- nrow(knownList); nu <- nrow(unknownList)
                #mk <- ncol(knownList)
                #if (nk>nu){NAList <- matrix(nrow=nk-nu,ncol=mk);
                #           unknownList <- rbind(as.matrix(unknownList),NAList)}
                #          else{NAList <- matrix(nrow=nu-nk,ncol=mk);
                #           knownList <- rbind(as.matrix(knownList),NAList)}
                #MList <- cbind(knownList,unknownList)
                #rownames(MList) <- NULL
                #M1 <- paste(MList[,2],MList[,4],MList[,5],MList[,6],sep="_")
                #M2 <- paste(MList[,8],MList[,10],MList[,11],MList[,12], sep="_")
                #write.csv(data.frame(Known=M1,unknown=M2), paste("niche_feature_",k,".csv",sep=""))
                #}



##making fake matrices:

lnormMat <- function(m, cv){
	sdlog <- sqrt(log(cv^2+1))
	meanlog <- log(m) - sdlog^2/2
	return(rlnorm(length(m), meanlog, sdlog))
}

gammaMat <- function(m, cv){
   shape = m*cv^2
   scale= 1/(cv^2)
   return(rgamma(length(m), shape=shape, scale=scale))
}

nbinomMat <- function(m, cv){
	size= m/(cv^2*m-1)
   return(rnbinom(n=length(m),size=size,mu=m))
}
</source-file>

