
  if(!file.exists(class.file.name)&
     file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
     file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
     file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file))) {
    
    print("Classifying leaf and wood points...")
    rf_predict(c2c.file=c2c.file, 
               # rf_model = rf_model,
               rf_model_path = rf_model_path,
               class.file.name=class.file.name)
    
    
  } 
  