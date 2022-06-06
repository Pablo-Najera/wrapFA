.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}

StartWelcomeMessage <- function(){
  paste(c("========================================================\n",
          "wrapFA Package",
          " [Version ", utils::packageDescription("wrapFA")$Version,
          "; ",utils::packageDescription("wrapFA")$Date, "]\n",
          "More information: https://github.com/Pablo-Najera/wrapFA\n", "\n",
          "Please cite the lavaan and MplusAutomation packages", "\n",
          "========================================================\n"),
        sep="")
}
