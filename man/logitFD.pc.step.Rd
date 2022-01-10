\name{logitFD.pc.step}
\alias{logitFD.pc.step}
\title{Functional Principal Component Logistic Regression by stepwise order}
\description{
Fit of the Functional Principal Component Logistic Regression model with Functional Principal Components and nonfunctional variables included in the model according their prediction ability by an authomatic stepwise selection method.
}
\usage{
logitFD.pc.step(Response, FDobj = list(), nonFDvars = NULL)
}
\arguments{
\item{Response}{Binary (numeric or character) vector of observations of the response variable.}

\item{FDobj}{List of functional objects from fda package with the curves of the predictor functional variables.}

\item{nonFDvars}{Matrix or data frame with the observations of non-functional variables.}
}
\value{
\item{glm.fit}{glm object of the fitted model. The object allows to use the summary() function for printing a summary of the fit, the anova() function to produce an analysis of variance table, and to extract useful features as coefficients, effects, fitted.values or residuals}
\item{Intercept}{Intercept estimated parameter}
\item{betalist}{List of functional objects (fdobj) of fda package with the estimated parameter functions. Each element of the list corresponds to the associated functional predictor located in the same position of FDobj. All methods of fd package can be used as the plot() function among others.}
\item{PC.variance}{List of data frames with explained variability of functional principal components of functional predictors. Each element of the list contains the acumulation variance matrix corresponding to the associated functional variable in the same position.}
\item{ROC.curve}{Object of the roc() function of the pROC package for prediction ability testing of the model. The object can be printed, plotted, or passed to many other functions}
}
\author{Carmen Lucia Reina <carmenlureina@gmail.com> and Manuel Escabias <escabias@ugr.es>}
\examples{
library(fda.usc)
data(aemet)
Temp<-aemet$temp$data
Prec<-exp(aemet$logprec$data)
StationsVars<-aemet$df[,c("ind","altitude","longitude","latitude")]
StationsVars$North<-c(1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,1,0,0,1,1,1,1,1,
0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,1,1,1)
TempMonth<-matrix(0,73,12)
PrecMonth<-matrix(0,73,12)
for (i in 1:nrow(TempMonth)){
  TempMonth[i,1]<-mean(Temp[i,1:31])
  PrecMonth[i,1]<-mean(Prec[i,1:31])
  TempMonth[i,2]<-mean(Temp[i,32:59])
  PrecMonth[i,2]<-mean(Prec[i,32:59])
  TempMonth[i,3]<-mean(Temp[i,60:90])
  PrecMonth[i,3]<-mean(Prec[i,60:90])
  TempMonth[i,4]<-mean(Temp[i,91:120])
  PrecMonth[i,4]<-mean(Prec[i,91:120])
  TempMonth[i,5]<-mean(Temp[i,121:151])
  PrecMonth[i,5]<-mean(Prec[i,121:151])
  TempMonth[i,6]<-mean(Temp[i,152:181])
  PrecMonth[i,6]<-mean(Prec[i,152:181])
  TempMonth[i,7]<-mean(Temp[i,182:212])
  PrecMonth[i,7]<-mean(Prec[i,182:212])
  TempMonth[i,8]<-mean(Temp[i,213:243])
  PrecMonth[i,8]<-mean(Prec[i,213:243])
  TempMonth[i,9]<-mean(Temp[i,244:273])
  PrecMonth[i,9]<-mean(Prec[i,244:273])
  TempMonth[i,10]<-mean(Temp[i,274:304])
  PrecMonth[i,10]<-mean(Prec[i,274:304])
  TempMonth[i,11]<-mean(Temp[i,305:334])
  PrecMonth[i,11]<-mean(Prec[i,305:334])
  TempMonth[i,12]<-mean(Temp[i,335:365])
  PrecMonth[i,12]<-mean(Prec[i,335:365])
  }
FourierBasis<-create.fourier.basis(rangeval = c(1,12),nbasis=7)
BsplineBasis<-create.bspline.basis(rangeval = c(1,12),nbasis=8)
TempMonth.fd<-Data2fd(argvals = c(1:12), y=t(TempMonth),basisobj = FourierBasis)
PrecMonth.fd<-Data2fd(argvals = c(1:12), y=t(PrecMonth),basisobj = BsplineBasis)
logitFD.pc.step(Response=StationsVars$North,FDobj=list(TempMonth.fd,PrecMonth.fd),
 nonFDvars = StationsVars[,c("altitude","longitude")])
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
\concept{Functional Data Analysis (FDA)}
\concept{Functional Principal Components Analysis (FPCA)}
\concept{logistic regression}
