#CREATE Subsampled point cloud and SURFACE MODELS
print(paste("Input gridded point cloud and estimate normals"))

print(paste("Processing", input_file))
termId<-run(paste(cloudcompare, # call Cloud Compare. The .exe file folder must be in the system PATH
                  "-SILENT",
                  "-C_EXPORT_FMT", "ASC", "-PREC", 6, #Set asc as export format
                  "-NO_TIMESTAMP",
                  "-COMPUTE_NORMALS",
                  "-O", input_file, #open the subsampled file
                  "-SAVE_CLOUDS",
                  sep = " "))

while (is.null(rstudioapi::terminalExitCode(termId))) {
  Sys.sleep(0.1)
}

result <- rstudioapi::terminalBuffer(termId)

# Delete the buffer and close the session in the IDE
rstudioapi::terminalKill(termId)

output_file = gsub(".ptx",".asc", input_file)



