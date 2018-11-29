getSampleDescription <- function(sample, print = FALSE){
  
  sampleMean = mean(sample)
  disp = var(sample)
  movedDisp = disp*(length(sample) - 1)/length(sample)
  median = median(sample)
  min = min(sample)
  max = max(sample)
  r = max - min
  variationCoef = 100*((movedDisp)^(1/2))/sampleMean
  oscCoef = 100*r/sampleMean
  
  if(print){
    print(paste("Sample mean: ", sampleMean))
    print(paste("Sample disp (fixed).: ", disp))
    print(paste("Moved disp.:", movedDisp))
    print(paste("Median: ", median))
    print(paste("Min: ", min))
    print(paste("Max: ", max))
    print(paste("Swing (R): ", r))
    print(paste("Variation coef.: ", variationCoef))
    print(paste("Oscillation coef.: ", oscCoef))
  }
  
  return(list("sampleMean" = sampleMean, "disp" = disp, "movedDisp" = movedDisp,
              "median" = median, "min" = min, "max" = max, "r" = r, "variationCoef" = variationCoef,
              "oscCoef" = oscCoef))
}

autoHist <- function(sample, print = FALSE){
  steps = (round(max(sample), digits = 1) - round(min(sample), digits = 1))*10 + 1

  par(mar = rep(2, 4))
  
  hist(testX, breaks = steps, freq = FALSE)
  lines(density(testX)) #плотность
  plot(ecdf(testX)) #эмпирическая функция распределения
}

confidenceIntervalEV <- function(sample, disp, confidenceLevel, print = FALSE){
  
  mean = mean(sample)
  stdiv = (disp)^(1/2)
  cGamma = qnorm((1 + confidenceLevel)/2) #Лаплас
  
  left = mean - (cGamma*stdiv)/((length(sample))^(1/2))
  right = mean + (cGamma*stdiv)/((length(sample))^(1/2))
  
  if(print){
    print("With known dispersion:")
    print(paste("Interval for expected value: (", left, "; ", right, ")", sep = ""))
    print(paste("With C = ", cGamma))
    print(paste("With stdiv = ", stdiv))
  }
  
  return(list("left" = left, "right" = right))
}

confidenceIntervalEVFull <- function(sample, confidenceLevel, print = FALSE){
  
  mean = mean(sample)
  stdiv = (var(sample))^(1/2)
  
  student = studentQuantile(confidenceLevel = confidenceLevel, n = length(sample) - 1, print = FALSE)
  
  left = mean - (student*stdiv)/((length(sample))^(1/2))
  right = mean + (student*stdiv)/((length(sample))^(1/2))
  
  if(print){
    print("With unknown dispersion:")
    print(paste("Interval for expected value: (", left, "; ", right, ")", sep = ""))
    print(paste("With t{n-1} = ", student))
    print(paste("With stdiv = ", stdiv))
  }
  
  return(list("left" = left, "right" = right))
  
}

confidenceIntervalDisp <- function(sample, confidenceLevel, print = FALSE){
  
  disp = var(sample)
  
  c = chisqCoefs(confidenceLevel = confidenceLevel, n = length(sample) - 1)
  c1 = c$c1
  c2 = c$c2
  
  left = (length(sample) - 1)*disp/c2
  right = (length(sample) - 1)*disp/c1
  
  if(print){
    print(paste("Confidence interval for dispersion:: (", left, "; ", right, ")", sep = ""))
    print(paste("With C1 = ", c1))
    print(paste("With C2 = ", c2))
  }
  
  return(list("left" = left, "right" = right))
}

normalDistributionCheck <- function(sample, alpha, print = FALSE){
  varSeries = table(sample)
  steps = (round(max(sample), digits = 1) - round(min(sample), digits = 1))*10 + 1
  p = numeric(steps)
  n = numeric(steps)
  
  stddiv = (var(sample)*(length(sample) - 1)/length(sample))^(1/2)
  mean = mean(sample)
  
  right = round(min(sample) + 0.049, digits = 1)
  val = (right - mean)/stddiv
  p[1] = pnorm(val)
  n[1] = sum(sample < right)
  
  for(i in 2:(steps - 1)){
    left = right
    right = right + 0.1
    p[i] = pnorm((right - mean)/stddiv) - pnorm((left - mean)/stddiv)
    n[i] = sum((sample >= left) & (sample < right))
  }
  
  p[6] = 1 - pnorm((right - mean)/stddiv)
  n[6] = sum(sample > right)
  
  
  
  length = length(sample)
  
  chisqn = 0
  for(i in 1:steps) chisqn = chisqn + ((n[i] - p[i]*length)^2)/(p[i]*length)
  
  crit = 0
  prob = 0
  while(prob < 1 - alpha){ #Перебираем точки, чтобы получить вероятность
    prob = pchisq(crit, df = steps - 3)
    crit = crit + 0.0001
  }
  
  if(print){
    print(n)
    print(p)
    print(paste("Crit stat = ", chisqn))
    print(paste("Critical dot = ", crit))
  }
  
  print("N agreement:")
  if(chisqn < crit) print("+")
  else print("-")
  
}

homogeneityCheck <- function(sample1, sample2, alpha, print = FALSE){
  
  #частоты
  steps1 = (round(max(sample1), digits = 1) - round(min(sample1), digits = 1))*10 + 1
  n1 = numeric(steps1)
  n2 = numeric(steps1)
  right = round(min(sample1) + 0.049, digits = 1)
  n1[1] = sum(sample1 < right)
  n2[1] = sum(sample2 < right)
  for(i in 2:(steps1 - 1)){
    left = right
    right = right + 0.1
    n1[i] = sum((sample1 >= left) & (sample1 < right))
    n2[i] = sum((sample2 >= left) & (sample2 < right))
  }
  n1[6] = sum(sample1 > right)
  n2[6] = sum(sample2 > right)
  
  n = n1 + n2
  
  length1 = length(sample1)
  length2 = length(sample2)
  sumLength = length1 + length2

  chisqn = 0
  
  for(i in 1:steps1){
    chisqn = chisqn + ((n1[i] - (n[i]*length1)/sumLength)^2)/(n[i]*length1) + ((n2[i] - (n[i]*length2)/sumLength)^2)/(n[i]*length2)
  }
  
  chisqn = chisqn*sumLength
  
  crit = 0
  prob = 0
  while(prob < 1 - alpha){
    prob = pchisq(crit, df = steps1 - 1)
    crit = crit + 0.0001
  }
  
  if(print){
    print(n1)
    print(n2)
    print(n)
    print(paste("Crit stat : ", chisqn))
    print(paste("Critical dot = ", crit))
  }
  
  print("Homogenity:")
  if(chisqn < crit) print("+")
  else print("-")
}

dispAndEVCheck <- function(sample, EV, disp, alpha, print = FALSE){
  
  stddiv = (var(sample))^(1/2)
  t = ((length(sample))^(1/2))*(mean(sample) - EV)/(stddiv)
  
  chisq = (var(sample) * (length(sample) - 1))/disp
  
  prob = 0
  student = 0
  
  #1 - (1 - conf)/2
  #conf ~ P
  #результат в виде t{alpha, k}, тест в R - односторонний
  
  while(prob < (1 - alpha/2)){
    prob = pt(student, df = length(sample) - 1, ncp = 0) #нужна обратная функция
    student = student + 0.0001
  }
  
  confidenceLevel = 1 - alpha
  
  c1 = 0
  prob = 0
  while(prob < (1 - confidenceLevel)/2){
    prob = pchisq(c1, df = length(sample) - 1)
    c1 = c1 + 0.0001
  }
  
  c2 = c1
  prob = 0
  while(prob < (1 + confidenceLevel)/2){
    prob = pchisq(c2, df = length(sample) - 1)
    c2 = c2 + 0.0001
  }
  
  if(print){
    print(paste("t = ", t))
    print(paste("Chisq = ", chisq))
    print(paste("Student", student))
    print(paste("With C1 = ", c1))
    print(paste("With C2 = ", c2))
  }
  
  print(paste("Is EV = ", EV))
  if(abs(t) > student) print("-")
  else print("+")
  
  print(paste("Is variance = ", disp))
  if((chisq > c1) && (chisq < c2)) print("+")
  else print("-")
}

meansEqualCheck <- function(x, y, alpha, print = FALSE){
  
  xmean = mean(x)
  ymean = mean(y)
  
  xdisp = var(x)
  ydisp = var(y)
  
  length = length(x)
  
  t = ((length/2)^(1/2))*(xmean - ymean)/((((length - 1)*xdisp + (length - 1)*ydisp)/(2*length - 2))^(1/2))
  
  prob = 0
  student = 0
  
  #1 - (1 - conf)/2
  #conf ~ P
  #результат в виде t{alpha, k}, тест в R - односторонний
  
  while(prob < (1 - alpha/2)){                  ###????????????????????
    prob = pt(student, df = 2*length(x) - 2, ncp = 0) #нужна обратная функция
    student = student + 0.0001
  }
  
  if(print){
    print(paste("Crit stat: ", t))
    print(paste("Crit dot: ", student))
  }
  
  print("Is EVx = EVy:")
  if(abs(t) < student) print("+")
  else print("-")
}

dispEqualCheck <- function(x, y, alpha, print = FALSE){
  
  xdisp = var(x)
  ydisp = var(y)
  
  f = xdisp/ydisp
  
  fs = fisherCoefs(confidenceLevel = 1 - alpha, n = length(sample) - 1, m = length(sample) - 1)
  f1 = fs$c1
  f2 = fs$c2
  
  if(print){
    print(paste("F = ", f))
    print(fs)
  }
  
  print("Is VarX = VarY:")
  if((f > f1) && (f < f2)) print("+")
  else print("-")
}