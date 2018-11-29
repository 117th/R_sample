chisqCoefs <- function(alpha = NULL, n, confidenceLevel = NULL, print = FALSE){
  
  if(!is.null(confidenceLevel)){
    
    c1 = 0
    prob = 0
    while(prob < (1 - confidenceLevel)/2){ #ХИ2.ОБР((1 - confidenceLevel)/2; level) - нижняя квантиль
      prob = pchisq(c1, n)
      c1 = c1 + 0.0001
    }
    
    c2 = c1
    prob = 0
    while(prob < (1 + confidenceLevel)/2){ #ХИ2.ОБР.ПХ((1 - confidenceLevel)/2) - верхняя
      prob = pchisq(c2, n)
      c2 = c2 + 0.0001
    }
    
    if(print){
      print(paste("With C1 = ", c1))
      print(paste("With C2 = ", c2))
    }
    return(list("c1" = c1, "c2" = c2))
  }
  else if(!is.null(alpha)){
    
    confidenceLevel = 1 - alpha
    
    c1 = 0
    prob = 0
    while(prob < alpha/2){
      prob = pchisq(c1, df = n)
      c1 = c1 + 0.0001
    }
    
    c2 = c1
    prob = 0
    while(prob < (2 - alpha)/2){
      prob = pchisq(c2, df = n)
      c2 = c2 + 0.0001
    }
    
    if(print){
      print(paste("With C1 = ", c1))
      print(paste("With C2 = ", c2))
    }
    return(list("c1" = c1, "c2" = c2))
  }
}

studentQuantile <- function(alpha = NULL, n, confidenceLevel = NULL, print = FALSE){
  
  if(!is.null(confidenceLevel)){
    
    prob = 0
    student = 0.1
    
    #1 - (1 - conf)/2
    #conf ~ P
    #результат в виде t{alpha, k}, тест в R - односторонний
    
    while(prob < 1 - (1 - confidenceLevel)/2){
      prob = pt(student, df = n, ncp = 0) #нужна обратная функция
      student = student + 0.0001
    }
    
    if(print){
      print(paste("With t = ", student))
    }
    return(student)
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