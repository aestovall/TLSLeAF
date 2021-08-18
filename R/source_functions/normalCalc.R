
  #Calculate normals for gridded TLS point cloud
  if(!file.exists(output_file)|
     overwrite) normalCalc(input_file)
  print("Done")
  