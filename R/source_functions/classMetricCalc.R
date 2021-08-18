
  
  #Classify wood and leaf from random forest classifier
  if(file.exists(angle.file.name)&
     !(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))|
     overwrite) classMetricCalc(angle.file.name, SS, scales)
  
  