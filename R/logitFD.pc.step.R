#' @name logitFD.pc.step
#' @rdname logitFD.pc.step
#'
#' @title Filtered Functional Principal Component Logistic Regression by stepwise order
#'
#' @description Fit of the Functional Principal Component Logistic Regression model with selected Functional Principal Components included in the model according their prediction ability.
#'
#' @param Response Binary (numeric or character) vector of observations of the response variable.
#' @param FDobj List of functional objects from fda package with the curves of the predictor functional variables.
#' @param nonFDvars Matrix or data frame with the observations of non-functional variables.
#'
#' @return  glm object of the fitted model, Intercept estimated parameter, list of functional objects with the estimated parameter functions, List of data frames with explained variability of functional principal components of functional predictors, and roc object of the pROC package for prediction ability testing of the model.
#'
#' @export
#' @import fda pROC expm
#' @importFrom stats glm princomp binomial step na.omit 
#' @importFrom graphics plot


logitFD.pc.step<-function(Response,FDobj=list(),nonFDvars=NULL){
  
  #escalar
  scalar<-list()
  for (i in 1:length(FDobj)){
    scalar[[i]]<-inprod(FDobj[[i]]$basis,FDobj[[i]]$basis)
  }
  
  #acp2
  pc<-list()
  for (i in 1:length(FDobj)){
    pc[[i]]<-pca.fd(FDobj[[i]],nharm=FDobj[[i]]$basis$nbasis)
  }
  
  #Response
  Response<-unlist(Response)
  Response<-as.vector(Response)
  
  #Step
  scores<-list()
  for (i in 1:length(FDobj)){
    scores[[i]]<-data.frame(pc[[i]]$scores)
    colnames(scores[[i]])<-c(paste(LETTERS[i],1:ncol(scores[[i]])))
    scores[[i]]<-data.frame(scores[[i]])
  }
  
  #Matriz variables funcionales
  if(!is.null(nonFDvars)){
    nonFDvars<-data.frame(nonFDvars)
    design<-data.frame(Response,scores,nonFDvars)
  }else {
    design<-data.frame(Response,scores)
  }
  
  fit.0<-glm(Response~1,data=design,family=binomial)
  fit.1<-glm(Response~.,data=design,family=binomial)
  fit.step<-step(fit.0,scope=list(lower=fit.0,upper=fit.1),direction="both",trace=0)
  design.fit<-design[names(fit.step$coefficients)[-1]]
  #managing one column dataframes...
  design.fit<-design.fit[,1:ncol(design.fit),drop=F]
  
  #nonFDcoefs
  if(!is.null(nonFDvars)){
    nonFDcoefs<-fit.step$coefficients[names(nonFDvars)]
    nonFDcoefs<-na.omit(nonFDcoefs)
  }else{
    nonFDcoefs<-NULL
  }
  
  #omitir variables no funcionales
  if(!is.null(nonFDvars)){
    design.fit<-design.fit[,!colnames(design.fit)==names(nonFDcoefs)]
  }else{
    design.fit<-design.fit
  }
  
  
  design.step<-list()
  for(i in 1:length(FDobj)){
    design.step[[i]]<-design.fit[grep(c(LETTERS[i]),names(design.fit))]
  }
  
  harmonics<-list()
  for(i in 1:length(FDobj)){
    harmonics[[i]]<-pc[[i]]$harmonics$coefs
    colnames(harmonics[[i]])<-c(paste(LETTERS[i],1:ncol(scores[[i]])))
    harmonics[[i]]<-data.frame(harmonics[[i]])
  }
  
  harmonics.step<-list()
  for(i in 1:length(FDobj)){
    harmonics.step[[i]]<-harmonics[[i]][names(design.step[[i]])]
  }
  
  #Funcion parametro
  beta<-list()
  beta.plot<-list()
  titles<-c()
  for(i in 1:length(FDobj)){
    titles[i]<-c(paste("Parameter function beta_",i))
  }
  for(i in 1:length(FDobj)){
    ifelse(ncol(harmonics.step[[i]])==1,beta[[i]]<-fd(as.matrix(harmonics.step[[i]])*fit.step$coefficients[names(design.step[[i]])],basisobj=FDobj[[i]]$basis),beta[[i]]<-fd(as.matrix(harmonics.step[[i]])%*%fit.step$coefficients[names(design.step[[i]])],basisobj=FDobj[[i]]$basis))
#    beta.plot[[i]]<-plot(beta[[i]],main=titles[i])
  }
  
  
  #Intercept
  #media coef.b?sicos
  gamma<-list()
  for(i in 1:length(FDobj)){
    gamma[[i]]<-inprod(mean.fd(FDobj[[i]]),beta[[i]])
  }
  Intercept<-as.numeric(fit.step$coefficients[1])-sum(unlist(gamma))
  
  #Tabla Variabilidad Explicada
  Comp<-list()
  for(i in 1:length(FDobj)){
    Comp[[i]]<-names(scores[[i]])
  }
  #Varianza explicada
  prop.var<-list()
  for(i in 1:length(FDobj)){
    prop.var[[i]]<-round((pc[[i]]$varprop)*100,3)
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
  
  
  #Curva ROC
  Predictor<-as.vector(fit.step$fitted.values)
  roc<-roc(Response,Predictor,percent=TRUE,auc=TRUE,quiet = T)
#  roc.plot<-plot.roc(roc,legacy.axes=TRUE,col="blue",main="Curva ROC",print.auc=TRUE)
  
  #Salida de resultados
  Results<-list(fit.step,nonFDcoefs,Intercept,beta,PC.var,roc)
  names(Results)<-c("glm.fit","nonFDcoefs","Intercept","betalist","PC.variance","ROC")
  return(invisible(Results))
  
  
  #Mostrar en pantalla
#  Results1<-list(beta.plot,roc.plot)
#  names(Results1)<-c("beta.plot","ROC.plot")
#  return(Results1)
  
}
