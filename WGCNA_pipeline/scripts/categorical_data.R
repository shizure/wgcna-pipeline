categorical_to_numeric = function(df) {
  # Importing the dataset
  dataset = df
  
  traits = dataset[2:49]
  
  # Encoding categorical data
  traits$Gender = factor(traits$Gender,
                         levels = c('female', 'male'),
                         labels = c(0, 1))
  
  traits$KDOQI = factor(traits$KDOQI,
                        levels = c('CKD 2 (Mild)', 'CKD 3 (Moderate', 'Normal kidney function', 'CKD 4 (Severe)'),
                        labels = c(1, 2, 3, 4))
  
  traits$BMI_WHO = factor(traits$BMI_WHO,
                          levels = c('Underweight', 'Normal', 'Overweight', 'Obese'),
                          labels = c(1, 2, 3, 4))
  
  traits$SmokerStatus = factor(traits$SmokerStatus,
                               levels = c('Never smoked', 'Ex-smoker', 'Current smoker'),
                               labels = c(1, 2, 3))
  
  traits$AlcoholUse = factor(traits$AlcoholUse,
                             levels = c('No', 'Yes'),
                             labels = c(0, 1))
  
  traits$DiabetesStatus = factor(traits$DiabetesStatus,
                                 levels = c('Control (no Diabetes Dx/Med)', 'Diabetes'),
                                 labels = c(0, 1))
  
  traits$Hypertension.selfreport = factor(traits$Hypertension.selfreport,
                                          levels = c('no', 'yes'),
                                          labels = c(0, 1))
  
  traits$Hypertension.selfreportdrug = factor(traits$Hypertension.selfreportdrug,
                                              levels = c('no', 'yes'),
                                              labels = c(0, 1))
  
  traits$Hypertension.composite = factor(traits$Hypertension.composite,
                                         levels = c('no', 'yes'),
                                         labels = c(0, 1))
  
  traits$Hypertension.drugs = factor(traits$Hypertension.drugs,
                                     levels = c('no', 'yes'),
                                     labels = c(0, 1))
  
  traits$Med.anticoagulants = factor(traits$Med.anticoagulants,
                                     levels = c('no', 'yes'),
                                     labels = c(0, 1))
  
  traits$Med.all.antiplatelet = factor(traits$Med.all.antiplatelet,
                                       levels = c('no', 'yes'),
                                       labels = c(0, 1))
  
  traits$Med.Statin.LLD = factor(traits$Med.Statin.LLD,
                                 levels = c('no', 'yes'),
                                 labels = c(0, 1))
  
  traits$Symptoms.5G = factor(traits$Symptoms.5G,
                              levels = c('Asymptomatic', 'Other', 'Ocular', 'Retinal infraction', 'TIA', 'Stroke'),
                              labels = c(0, 1, 2, 3, 4, 5))
  
  traits$AsymptSympt = factor(traits$AsymptSympt,
                              levels = c('Asymptomatic', 'Ocular and others', 'Symptomatic'),
                              labels = c(0, 1, 2))
  
  traits$AsymptSympt2G = factor(traits$AsymptSympt2G,
                                levels = c('Asymptomatic', 'Symptomatic'),
                                labels = c(0, 1))
  
  traits$CAD_history = factor(traits$CAD_history,
                              levels = c('No history CAD', 'History CAD'),
                              labels = c(0, 1))
  
  traits$PAOD = factor(traits$PAOD,
                       levels = c('no', 'yes'),
                       labels = c(0, 1))
  
  traits$EP_composite = factor(traits$EP_composite,
                               levels = c('No composite endpoints', 'Composite endpoints'),
                               labels = c(0, 1))
  
  traits$EP_major = factor(traits$EP_major,
                           levels = c('No major events (endpoints)', 'Major events (endpoints)'),
                           labels = c(0, 1))
  
  traits$Macrophages.bin = factor(traits$Macrophages.bin,
                                  levels = c('no/minor', 'moderate/heavy'),
                                  labels = c(0, 1))
  
  traits$SMC.bin = factor(traits$SMC.bin,
                          levels = c('no/minor', 'moderate/heavy'),
                          labels = c(0, 1))
  
  traits$IPH.bin = factor(traits$IPH.bin,
                          levels = c('no', 'yes'),
                          labels = c(0, 1))
  
  traits$Calc.bin = factor(traits$Calc.bin,
                           levels = c('no/minor', 'moderate/heavy'),
                           labels = c(0, 1))
  
  traits$Collagen.bin = factor(traits$Collagen.bin,
                               levels = c('no/minor', 'moderate/heavy'),
                               labels = c(0, 1))
  
  traits$Fat.bin_10 = factor(traits$Fat.bin_10,
                             levels = c('<10%', '>10%'),
                             labels = c(0, 1))
  
  traits$Fat.bin_40 = factor(traits$Fat.bin_40,
                             levels = c('<40%', '>40%'),
                             labels = c(0, 1))
  
  traits$OverallPlaquePhenotype = factor(traits$OverallPlaquePhenotype,
                                         levels = c('fibrous', 'atheromatous', 'fibroatheromatous'),
                                         labels = c(0, 1, 2))
  
  #######
  
  #Convert all columns to numeric values
  for(i in names(traits)){
    if(!(all(is.na(traits[[i]])))){
      traits[[i]] = as.numeric(as.character(traits[[i]]))
    }
  }
  return(cbind(dataset[1], traits))
}
