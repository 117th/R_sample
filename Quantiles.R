chisqCoefs <- function(alpha = NULL, n, confidenceLevel = NULL, print = FALSE){
  
  c1 = qchisq((1 - confidenceLevel)/2, n - 1, lower.tail = TRUE)
  c2 = qchisq((1 + confidenceLevel)/2, n - 1, lower.tail = TRUE)
    
  if(print){
    print(paste("With qc1 = ", qchisq((1 - confidenceLevel)/2, n - 1, lower.tail = TRUE)))
    print(paste("With qc2 = ", qchisq((1 + confidenceLevel)/2, n - 1, lower.tail = TRUE)))
  }
  
  return(list("c1" = c1, "c2" = c2))
  
}

studentQuantile <- function(alpha = NULL, n, confidenceLevel = NULL, print = FALSE){
  
  if(!is.null(confidenceLevel)){
    if(print){
      print(paste("qt =", qt((1 - confidenceLevel)/2, n, lower.tail = FALSE)))
    }
    return(qt((1 - confidenceLevel)/2, n, lower.tail = FALSE))
  }
  else if(!is.null(alpha)){
    
    prob = 0
    student = 0.1
    
    while(prob < 1 - alpha/2){
      prob = pt(student, df = n, ncp = 0) #нужна обратная функция
      student = student + 0.0001
    }
    
    if(print){
      print(paste("With t = ", student))
    }
    return(student)
  }
}

fisherCoefs <- function(alpha = NULL, n, m, confidenceLevel = NULL, print = FALSE){
  
  if(!is.null(confidenceLevel)){
    c1 = 0
    prob = 0
    while(prob < (1 - confidenceLevel)/2){ #ХИ2.ОБР((1 - confidenceLevel)/2; level) - нижняя квантиль
      prob = pf(c1, n, m, ncp = 0)
      c1 = c1 + 0.0001
    }
    
    c2 = c1
    prob = 0
    while(prob < (1 + confidenceLevel)/2){ #ХИ2.ОБР.ПХ((1 - confidenceLevel)/2) - верхняя
      prob = pf(c2, n, m, ncp = 0)
      c2 = c2 + 0.0001
    }
    
    if(print){
      print(paste("With C1 = ", c1))
      print(paste("With C2 = ", c2))
    }
    return(list("c1" = c1, "c2" = c2))
  }
  
}