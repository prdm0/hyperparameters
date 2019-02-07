config_cpp <- function(dir = NULL){
   
   ifelse(is.null(dir), base::setwd(base::getwd()), 
          base::setwd(dir))
   
   if(!("fast.cpp" %in% base::dir())){
      stop("The file fast.cpp is not found in the directory: ",
           base::getwd())
   }
   
   depend <- function(...){
      return("Rcpp" %in% utils::installed.packages(...)[ , "Package"])
   }
   
   if(depend()){
      Rcpp::sourceCpp("fast.cpp")
      return(message("All ready. Setup completed!"))
   }else{
      message("==> The Rcpp package will need to be installed.")
      install <- function(...){
         tryCatch(expr = utils::install.packages(...),
                  warning = function(w) NA)
      } 
      
      pkg <- NULL
      pkg <- install("Rcpp") 
      
      if(!(is.null(pkg))) stop(" ==> Check your internet connection!")
      
      Rcpp::sourceCpp("fast.cpp")
      return(message("All ready. Setup completed!"))
   }
}
