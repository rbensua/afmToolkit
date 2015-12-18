afmBaselineCorrection <- function(afmdata, ZPointApp = NULL, ZPointRet = NULL){
  # First determine how many segments are in the curve
  if ("CP" %in% names(afmdata) & is.null(ZPointApp)){
    Zapp <- subset(afmdata$data,Segment = "approach")$Z  
    ZPointApp <- 0.3*max(Zapp)+0.7*afmdata$CP[["CP"]]  
  }
  if (is.null(ZPointRet)){
    ZPointRet = ZPointApp 
  }

    N <- nlevels(afmdata$data$Segment)
  data.approach <- subset(afmdata$data, Segment == "approach" & 
                          abs(Z) > abs(ZPointApp), 
                          select = c("Z","Force"))
  fit.approach <- lm(Force ~ Z, data = data.approach)
  F.corrected.approach <- subset(afmdata$data, Segment == "approach")$Force - 
    predict(fit.approach,data.frame(Z = subset(afmdata$data,
                                               Segment == "approach")$Z))
  
  if (N == 1){
    # If N = 1 there is only the approach  
#    ForceCorrected <- list(F.corrected.approach = F.corrected.approach)
    afmdata$data$ForceCorrected <- c(F.corrected.approach)
  }else if (N == 2){
    if ("CP" %in% names(afmdata)){
      Zret <- subset(afmdata$data,Segment = "retract")$Z
      ZPointRet <- 0.6*max(Zret)+0.4*afmdata$CP[["CP"]]  
    }
    # If N = 2 there are approach and retract
    
    data.retract <- subset(afmdata$data, Segment == "retract" & 
                            abs(Z) > abs(ZPointRet), 
                            select = c("Z","Force"))
    fit.retract <- lm(Force ~ Z, data = data.retract)
    F.corrected.retract <- subset(afmdata$data, Segment == "retract")$Force - 
      predict(fit.retract,data.frame(Z = subset(afmdata$data, Segment == "retract")$Z))

    #    ForceCorrected <- list(F.corrected.approach = F.corrected.approach,
#                           F.corrected.retract = F.corrected.retract)
    afmdata$data$ForceCorrected <- c(F.corrected.approach,
                                     F.corrected.retract)
  }else{
    # If N = 3 there are approach, contact and retract.
    
    data.retract <- subset(afmdata$data, Segment == "retract" & 
                             abs(Z) > abs(ZPointRet), 
                           select = c("Z","Force"))
    
    F.corrected.contact <-subset(afmdata$data, Segment == "contact")$Force - 
      predict(fit.approach,data.frame(Z = subset(afmdata$data,
                                                 Segment == "contact")$Z))
    
    fit.retract <- lm(Force ~ Z, data = data.retract)
    F.corrected.retract <- subset(afmdata$data, Segment == "retract")$Force - 
      predict(fit.retract,data.frame(Z = subset(afmdata$data, Segment == "retract")$Z))
#     ForceCorrected <- list(F.corrected.approach = F.corrected.approach,
#                            F.corrected.contact = F.corrected.contact,
#                            F.corrected.retract = F.corrected.retract)
    afmdata$data$ForceCorrected <- c(F.corrected.approach,
                                     F.corrected.contact,
                                     F.corrected.retract)
  }
  
  return(afmdata(afmdata))
}