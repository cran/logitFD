#' @name logitFD.fpc
#' @rdname logitFD.fpc
#'
#' @title Filtered Functional Principal Component Logistic Regression by explained variability order
#'
#' @description Fit of the Filtered Functional Principal Component Logistic Regression model with selected Functional Principal Components included in the model according their explained variability.
#'
#' @param Response Binary (numeric or character) vector of observations of the response variable.
#' @param FDobj List of functional objects from fda package with the curves of the predictor functional variables.
#' @param ncomp Numeric vector with the number of components to be considered for each functional predictor. The vector has equal lenght than FDobj.
#' @param nonFDvars Matrix or data frame with the observations of non-functional variables.
#'
#' @return  glm object of the fitted model, Intercept estimated parameter, list of functional objects with the estimated parameter functions, List of data frames with explained variability of functional principal components of functional predictors, and roc object of the pROC package for prediction ability testing of the model.
#'
#'
#' @export
#' @import fda expm pROC
#' @importFrom stats glm princomp binomial na.omit 
#' @importFrom graphics plot

logitFD.fpc<-function(Response,FDobj=list(),ncomp=c(),nonFDvars=NULL){
  
  #escalar
  scalar<-list()
  for (i in 1:length(FDobj)){
    scalar[[i]]<-inprod(FDobj[[i]]$basis,FDobj[[i]]$basis)
  }
  
  #acp1
  fpc<-list()
  for (i in 1:length(FDobj)){
    fpc[[i]]<-pca.fd(fd(sqrtm(scalar[[i]])%*%FDobj[[i]]$coefs,basisobj=FDobj[[i]]$basis),nharm=FDobj[[i]]$basis$nbasis)
  }
  
  #Response
  Response<-unlist(Response)
  Response<-as.vector(Response)
  
  for(i in 1:length(FDobj)){
    if(ncomp[i]>FDobj[[i]]$basis$nbasis){
      ncomp[i]<-c(FDobj[[i]]$basis$nbasis)
      warning("the number of components can?t exceeed the number of basis")
    }else{
      ncomp<-ncomp
    }
  }
  
  #Ajuste
  scores<-list()
  scores.var<-list()
  for (i in 1:length(FDobj)){
    scores[[i]]<-fpc[[i]]$scores
    colnames(scores[[i]])<-c(paste(LETTERS[i],1:ncol(scores[[i]])))
    scores[[i]]<-data.frame(scores[[i]])
    scores.var[[i]]<-scores[[i]][,1:ncomp[i],drop=F]
  }
  
  #Matriz variables funcionales
  if(!is.null(nonFDvars)){
    nonFDvars<-data.frame(nonFDvars)
    design<-data.frame(Response,scores.var,nonFDvars)
  }else {
    design<-data.frame(Response,scores.var)
  }
  
  #Ajuste
  fit<-glm(design,family=binomial)
  design.fit<-design[names(fit$coefficients)[-1]]
  
  
  #nonFDcoefs
  if(!is.null(nonFDvars)){
    nonFDcoefs<-fit$coefficients[names(nonFDvars)]
  }else{
    nonFDcoefs<-NULL
  }
  
  #quedan eliminadas las variables no funcionales
  design.var<-list()
  for(i in 1:length(FDobj)){
    design.var[[i]]<-design[names(scores.var[[i]])]
  }
  
  harmonics<-list()
  for(i in 1:length(FDobj)){
    harmonics[[i]]<-fpc[[i]]$harmonics$coefs
    colnames(harmonics[[i]])<-c(paste(LETTERS[i],1:ncol(scores[[i]])))
    harmonics[[i]]<-data.frame(harmonics[[i]])
  }
  
  harmonics.var<-list()
  for(i in 1:length(FDobj)){
    harmonics.var[[i]]<-harmonics[[i]][names(design.var[[i]])]
  }
  
  #Funcion parametro 
  beta<-list()
  beta.plot<-list()
  titles<-c()
  for(i in 1:length(FDobj)){
    titles[i]<-c(paste("Parameter function beta_",i))
  }
  for(i in 1:length(FDobj)){
    ifelse(ncol(harmonics.var[[i]])==1,beta[[i]]<-fd(sqrtm(scalar[[i]])%*%as.matrix(harmonics.var[[i]])*fit$coefficients[names(design.var[[i]])],basisobj=FDobj[[i]]$basis),beta[[i]]<-fd(sqrtm(scalar[[i]])%*%as.matrix(harmonics.var[[i]])%*%fit$coefficients[names(design.var[[i]])],basisobj=FDobj[[i]]$basis))
#    beta.plot[[i]]<-plot(beta[[i]],main=titles[i])
  }
  
  #Intercept
  #media coef.b?sicos
  gamma<-list()
  for(i in 1:length(FDobj)){
    gamma[[i]]<-inprod(mean.fd(FDobj[[i]]),beta[[i]])
  }
  Intercept<-as.numeric(fit$coefficients[1])-sum(unlist(gamma))
  
  #Tabla Variabilidad Explicada
  Comp<-list()
  for(i in 1:length(FDobj)){
    Comp[[i]]<-names(scores[[i]])
  }
  #Varianza explicada
  prop.var<-list()
  for(i in 1:length(FDobj)){
    prop.var[[i]]<-round((fpc[[i]]$varprop)*100,3)
  }
  #Varianza acumulada
  cum.prop.var<-list()
  for(i in 1:length(FDobj)){
    cum.prop.var[[i]]<-cumsum(prop.var[[i]])
  }
  #Tabla resumen
  PC.var<-list()
  for(i in 1:length(FDobj)){
    PC.var[[i]]<-data.frame(Comp[[i]],prop.var[[i]],cum.prop.var[[i]])
    names(PC.var[[i]])<-c("Comp.","% Prop.Var","% Cum.Prop.Var")
  }
  
  #Curva ROC**
  Predictor<-as.vector(fit$fitted.values)
  roc<-roc(Response,Predictor,percent=TRUE,auc=TRUE,quiet = T)
#  roc.plot<-plot.roc(roc,legacy.axes=TRUE,col="blue",main="Curva ROC",print.auc=TRUE)
  
  #Resultados de salida
  Results<-list(fit,nonFDcoefs,Intercept,beta,PC.var,roc)
  names(Results)<-c("glm.fit","nonFDcoefs","Intercept","betalist","PC.variance","ROC")
  return(invisible(Results))
  
  #Mostrar en pantalla
#  Results1<-list(beta.plot,roc.plot)
#  names(Results1)<-c("beta.plot","ROC.plot")
#  return(Results1)
}
