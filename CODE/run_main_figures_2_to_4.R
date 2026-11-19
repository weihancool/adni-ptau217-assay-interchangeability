.revision_root <- Sys.getenv("STUDY_B_ROOT", "")
if (!nzchar(.revision_root)) {
  .revision_frames <- lapply(sys.frames(), function(x) x$ofile)
  .revision_files <- unlist(.revision_frames)
  if (!length(.revision_files)) {
    .revision_files <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
  }
  if (length(.revision_files)) {
    .revision_root <- sub("[/\\\\]CODE[/\\\\].*$", "", tail(.revision_files,1), useBytes=TRUE)
    .revision_root <- normalizePath(.revision_root,winslash="/")
  } else {
    .revision_root <- normalizePath(getwd(),winslash="/")
    for (.revision_i in 1:6) {
      if (file.exists(file.path(.revision_root,"CODE","42_final_display.R"))) break
      .revision_root <- dirname(.revision_root)
    }
  }
}
stopifnot(file.exists(file.path(.revision_root,"CODE","42_final_display.R")))
Sys.setenv(STUDY_B_ROOT=.revision_root)
source(file.path(.revision_root,"CODE","42_final_display.R"),
       local=new.env(parent=globalenv()),encoding="UTF-8")
