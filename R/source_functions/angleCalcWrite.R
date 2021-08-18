
  #Calculate scattering angle and leaf angle
  if((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
    rstudioapi::terminalBusy(rstudioapi::terminalList())])==0)&
    (!file.exists(angle.file.name)|
     overwrite)){
    # remove(dat)
    # gc()
    
    print("Calculating surface angles...")
    angleCalcWrite(output_file, center, scatterLim, angle.file.name)
    
  }
  