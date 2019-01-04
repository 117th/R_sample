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
  steps = 6

  par(mar = rep(2, 4))
  
  plot(hist(sample, breaks = 6, freq = FALSE))
  lines(density(sample), col = "red") #плотность
  plot(ecdf(sample)) #эмпирическая функция распределения
}

confidenceIntervalEV <- function(sample, disp, confidenceLevel, print = FALSE){
  
  mean = mean(sample)
  stdiv = (disp)^(1/2)
  cGamma = qnorm((1 + confidenceLevel)/2) #Лаплас
  
  left = mean - (cGamma*stdiv)/((length(sample))^(1/2))
  right = mean + (cGamma*stdiv)/((length(sample))^(1/2))
  
  if(print){
    print("With known dispersion:")
    print(paste("Interval for expected value: (", left, ";  ", right, ")", sep = ""))
    print(paste("With C = ", cGamma))
    print(paste("With stdiv = ", stdiv))
  }
  
  return(list("interval" = paste("(", round(left, digits = 3), "; ", round(right, digits = 3), ")", sep =""), "cGamma" = cGamma))
}

confidenceIntervalEVFull <- function(sample, confidenceLevel, print = FALSE){
  
  mean = mean(sample)
  stdiv = (var(sample))^(1/2)
  
  student = qt((1 - confidenceLevel)/2, df = length(sample) - 1, lower.tail = FALSE)
  
  left = mean - (student*stdiv)/((length(sample))^(1/2))
  right = mean + (student*stdiv)/((length(sample))^(1/2))
  
  if(print){
    print("With unknown dispersion:")
    print(paste("Interval for expected value: (", left, "; ", right, ")", sep = ""))
    print(paste("With t{n-1} = ", student))
    print(paste("With stdiv = ", stdiv))
  }
  
  return(list("interval" = paste("(", round(left, digits = 3), "; ", round(right, digits = 3), ")", sep =""), 
              "student" = student))
  
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
  
  return(list("interval" = paste("(", round(left, digits = 3), "; ", round(right, digits = 3), ")", sep =""), 
              "c1" = round(c1, digits = 3), "c2" = round(c2, digits = 3)))
}

normalDistributionCheck <- function(sample, alpha, print = FALSE){
  varSeries = table(sample)
  steps = (round(max(sample), digits = 1) - round(min(sample), digits = 1))*10 + 1
  p = numeric(steps)
  n = numeric(steps)
  names = c()
  
  stddiv = (var(sample)*(length(sample) - 1)/length(sample))^(1/2)
  mean = mean(sample)
  
  right = round(min(sample) + 0.049, digits = 1)
  val = (right - mean)/stddiv
  p[1] = pnorm(val)
  n[1] = sum(sample < right)
  names[1] = paste("(-inf;", right, ")")
  
  for(i in 2:(steps - 1)){
    left = right
    right = right + 0.1
    p[i] = pnorm((right - mean)/stddiv) - pnorm((left - mean)/stddiv)
    n[i] = sum((sample >= left) & (sample < right))
    names[i] = paste("[", left, ";", right, ")")
  }
  
  p[6] = 1 - pnorm((right - mean)/stddiv)
  n[6] = sum(sample > right)
  
  length = length(sample)
  
  chisqn = 0
  for(i in 1:steps) chisqn = chisqn + ((n[i] - p[i]*length)^2)/(p[i]*length)
 
  crit = qchisq(1 - alpha,  steps - 3, lower.tail = TRUE)
  
  if(print){
    print(names)
    print(n)
    print(p)
    print(paste("With delta =", steps))
    print(paste("Crit stat = ", chisqn))
    print(paste("Critical dot = ", crit))
  }
  
  print("N agreement:")
  if(chisqn < crit) print("+")
  else print("-")
  
}

homogeneityCheck <- function(sample1, sample2, alpha, print = FALSE){
  
  #частоты
  steps1 = 6
  incr = (round(max(sample1), digits = 1) - round(min(sample1), digits = 1))/6
  n1 = numeric(steps1)
  n2 = numeric(steps1)
  names = c()
  
  right = round(min(sample1) + 0.049, digits = 1)
  n1[1] = sum(sample1 < right)
  n2[1] = sum(sample2 < right)
  names[1] = paste("(-inf;", right, ")")
  for(i in 2:(steps1 - 1)){
    left = right
    right = right + incr
    n1[i] = sum((sample1 >= left) & (sample1 < right))
    n2[i] = sum((sample2 >= left) & (sample2 < right))
    names[i] = paste("[", left, ";", right, ")")
  }
  n1[6] = sum(sample1 > right)
  n2[6] = sum(sample2 > right)
  names[6] = paste("[", right, "; +inf)")
  
  n = n1 + n2
  
  length1 = length(sample1)
  length2 = length(sample2)
  sumLength = length1 + length2

  chisqn = 0
  
  for(i in 1:steps1){
    chisqn = chisqn + ((n1[i] - (n[i]*length1)/sumLength)^2)/(n[i]*length1) + ((n2[i] - (n[i]*length2)/sumLength)^2)/(n[i]*length2)
  }
  
  chisqn = chisqn*sumLength
  
  crit = qchisq(1 - alpha,  steps1 - 1, lower.tail = TRUE)
  
  df = matrix(n1, ncol = steps1)
  df = rbind(df, n2)
  df = rbind(df, n)
  
  df = as.data.frame(df)
  
  if(print){
    print(names)
    print(n1)
    print(n2)
    print(n)
    print(paste("Crit stat : ", chisqn))
    print(paste("Critical dot = ", crit))
  }
  
  print("Homogenity:")
  if(chisqn < crit) print("+")
  else print("-")
  
  if(chisqn < crit) res = "H0 is ok"
  else res = "H0 is not ok"
  
  return(list("df" = df, "res" = res, 
              names = paste(names, collapse = ", "), 
              "stat" = round(chisqn, digits = 2), "crit" = round(crit, digits = 2)))
}

dispAndEVCheck <- function(sample, EV, disp, alpha, print = FALSE){
  
  stddiv = (var(sample))^(1/2)
  t = ((length(sample))^(1/2))*(mean(sample) - EV)/(stddiv)
  
  chisq = (var(sample) * (length(sample) - 1))/disp
  
  student = qt(alpha/2, length(sample) - 1, lower.tail = FALSE)
  
  c2 = qchisq(1 - alpha/2, df = length(sample) - 1)
  c1 = qchisq(alpha/2, df = length(sample) - 1)
  
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
  
  if(abs(t) > student) resEV = "Not equal (|t| > crit)"
  else resEV = "Equal (|t| < crit)"
  
  if((chisq > c1) && (chisq < c2)) resDisp = "Equal (stat in interval (c1; c2))"
  else resDisp = "Not equal (stat not in interval (c1; c2))"
  
  return(list("stat" = t, "crit" = student,
              "chisq" = chisq, "resDisp" = resDisp,
              "resEV" = resEV, "c1" = round(c1, digits = 3), 
              "c2" = round(c2, digits = 3)))
}

meansEqualCheck <- function(x, y, alpha, print = FALSE){
  
  xmean = mean(x)
  ymean = mean(y)
  
  xdisp = var(x)
  ydisp = var(y)
  
  length = length(x)
  
  t = ((length/2)^(1/2))*(xmean - ymean)/((((length - 1)*xdisp + (length - 1)*ydisp)/(2*length - 2))^(1/2))
  
  student = qt(1 - alpha/2, df = 2*length(x) - 2, ncp = 0)
  
  if(print){
    print(paste("Crit stat: ", t))
    print(paste("Crit dot: ", student))
  }
  
  print("Is EVx = EVy:")
  if(abs(t) < student) print("+")
  else print("-")
  
  if(abs(t) < student) res = "Mean1 = Mean2"
  else res = "Mean1 != Mean2"
  
  return(list("stat" = t, "crit" = student, "res" = res))
}

dispEqualCheck <- function(x, y, alpha, print = FALSE){
  
  xdisp = var(x)
  ydisp = var(y)
  
  f = xdisp/ydisp
  
  f1 = qf(alpha/2, df1 = length(sample) - 1, df2 = length(sample) - 1)
  f2 = qf((2 - alpha)/2, df1 = length(sample) - 1, df2 = length(sample) - 1)
  
  if(print){
    print(paste("F = ", f))
    print(paste("f1 =", f1))
    print(paste("f2 =", f2))
  }
  
  print("Is VarX = VarY:")
  if((f > f1) && (f < f2)) print("+")
  else print("-")
  
  if((f > f1) && (f < f2)) res = "Disp1 = Disp2"
  else res = "Disp1 != Disp2"
  
  return(list("f" = f, "f1" = f1, "f2" = f2, "res" = res))
}

normalDistributionCheck2 <- function(sample, alpha, print = FALSE){
  varSeries = table(sample)
  steps = 6
  incr = (round(max(sample), digits = 1) - round(min(sample), digits = 1))/6
  p = numeric(steps)
  n = numeric(steps)
  names = c()
  
  stddiv = (var(sample)*(length(sample) - 1)/length(sample))^(1/2)
  mean = mean(sample)
  
  print(paste("Stddiv =", stddiv))
  print(paste("Mean =", mean))
  
  right = round(min(sample) + 0.049, digits = 1)
  val = (right - mean)/stddiv
  p[1] = pnorm(val)
  n[1] = sum(sample < right)
  names[1] = paste("(-inf;", right, ")")
  
  for(i in 2:(steps - 1)){
    left = right
    right = right + incr
    p[i] = pnorm((right - mean)/stddiv) - pnorm((left - mean)/stddiv)
    n[i] = sum((sample >= left) & (sample < right))
    names[i] = paste("[", left, ";", right, ")")
  }
  
  p[6] = 1 - pnorm((right - mean)/stddiv)
  n[6] = sum(sample > right)
  names[6] = paste("[", right, "; +inf)")
  
  length = length(sample)
  
  chisqn = 0
  for(i in 1:steps) chisqn = chisqn + ((n[i] - p[i]*length)^2)/(p[i]*length)
  
  crit = qchisq(1 - alpha,  steps - 3, lower.tail = TRUE)
  
  if(print){
    print(names)
    print(n)
    print(p)
    print(paste("With delta =", incr))
    print(paste("Crit stat = ", chisqn))
    print(paste("Critical dot = ", crit))
  }
  
  df = matrix(n, ncol = steps)
  df = rbind(df, p)
  
  df = as.data.frame(df)
  #colnames(df) = make.names(c(1:steps))
  
  print("N agreement:")
  if(chisqn < crit) print("+")
  else print("-")
  
  if(chisqn < crit) res = "H0 is ok"
  else res = "H0 is not ok"
  
  return(list("df" = df, "res" = res, 
              names = paste(names, collapse = ", "), 
              "stat" = round(chisqn, digits = 2), "crit" = round(crit, digits = 2)))
}

corellationCheck <- function(sample, sample2, confidenceLevel, alpha, print = FALSE){
  
  meanX = mean(sample)
  meanY = mean(sample2)
  
  dispX = var(sample)*(length(sample) - 1)/length(sample)
  dispY = var(sample2)*(length(sample2) - 1)/length(sample2)
  
  length = length(sample)
  
  correlation = 0
  for(i in (1:length)) correlation = correlation + sample[i]*sample2[i]
  
  correlation = correlation  - length*meanY*meanX
  correlation = correlation/(length*(dispY*dispX)^(1/2))
  
  if(print){
    print(paste("rXY =", correlation))
  }
  
  lambda1 = 0.5*log((1 + correlation)/(1 - correlation)) - qnorm((1 + confidenceLevel)/2)/(length - 3)^(1/2)
  lambda2 = 0.5*log((1 + correlation)/(1 - correlation)) + qnorm((1 + confidenceLevel)/2)/(length - 3)^(1/2)

  left = tanh(lambda1)
  right = tanh(lambda2)
  
  if(print){
    print(paste("Interval for rXY: (", left, "; ", right, ")", sep = ""))
  }
  
  tStat = correlation*((length - 2)/(1 - correlation^2))^(1/2)
  tCrit = qt(1 - alpha/2, length - 2)
  
  if(print){
    print(paste("t stat =", tStat))
    print(paste("t crit =", tCrit))
  }
  
  print("Is rXY = 0")
  if(abs(tStat) < tCrit) print("+")
  else("-")
  
  if(abs(tStat) < tCrit) res = "r = 0"
  else res = "r != 0"
  
  return(list("corrCoef" = correlation, "stat" = round(tStat, digits = 3),
              "crit" = round(tCrit, digits = 3),
              "res" = res,
              "interval" = paste("(", round(left, digits = 3), "; ", round(right, digits = 3), ")", sep =""),
              "cGamma" = qnorm((1 + confidenceLevel)/2)))
}

probAnalysis <- function(sample, sample2, p0, a, k, sigma, confidenceLevel = 0.95, alpha = 0.05, print = FALSE){
  
  upperLim = a + k*sigma
  
  length1 = length(sample)
  length2 = length(sample2)
  
  pX = (sum(sample < upperLim))/length1
  
  cGamma = qnorm((1 + confidenceLevel)/2) #Лаплас
  
  left = pX - qnorm((1 + confidenceLevel)/2)*(pX*(1 - pX)/length1)^(1/2)
  right = pX + qnorm((1 + confidenceLevel)/2)*(pX*(1 - pX)/length1)^(1/2)
  
  if(print){
    print(paste("Interval for pX: (", left, "; ", right, ")", sep = ""))
    print(paste("pX =", pX))
  }
  
  pY = (sum(sample2 < upperLim))/length2
  
  left = pY - qnorm((1 + confidenceLevel)/2)*(pY*(1 - pY)/length2)^(1/2)
  right = pY + qnorm((1 + confidenceLevel)/2)*(pY*(1 - pY)/length2)^(1/2)
  
  if(print){
    print(paste("Interval for pY: (", left, "; ", right, ")", sep = ""))
    print(paste("pY =", pY))
  }
  
  pXInterval = paste("(", round(left, digits = 3), "; ", round(right, digits = 3), ")", sep ="")
  
  difference = pX - pY
  
  left = difference - qnorm((1 + confidenceLevel)/2)*((pX*(1 - pX) + pY*(1 - pY))/length1)^(1/2)
  right = difference + qnorm((1 + confidenceLevel)/2)*((pX*(1 - pX) + pY*(1 - pY))/length1)^(1/2)
  
  if(print){
    print(paste("Interval for pX - pY: (", left, "; ", right, ")", sep = ""))
    print(paste("pX - pY =", difference))
  }
  
  pXYInterval = paste("(", round(left, digits = 3), "; ", round(right, digits = 3), ")", sep ="")
  
  uStat = (pX - p0)/((p0*(1 - p0)/length1)^(1/2))
  uCrit = qnorm(1 - alpha/2)

  if(print){
    print("")
    print("Check pX ~ p0")
    print(paste("u stat =", uStat))
    print(paste("u crit =", uCrit))
  }
  
  if(abs(uStat) > uCrit) print("pX =!= p0")
  else print("pX == p0")
  
  if(abs(uStat) > uCrit) res0 = "pX =!= p0"
  else res0 = "pX = p0"
  
  uStatO = uStat
  uStat = (difference*(2*length1)^(1/2))/(((pX + pY) * (2 - pX - pY))^(1/2))
  
  if(print){
    print("")
    print("Check pX ~ pY")
    print(paste("u stat =", uStat))
    print(paste("u crit =", uCrit))
  }
  
  if(abs(uStat) > uCrit) print("pX =!= pY")
  else print("pX == pY")
  
  if(abs(uStat) > uCrit) res = "pX =!= pY"
  else res = "pX = pY"
  
  return(list("pX" = pX, "pY" = pY, "pXInterval" = pXInterval, "cGamma" = cGamma,
              "pXYInterval" = pXYInterval, "uStat" = uStat, "uStat0" = uStatO, "uCrit" = uCrit, "res0" = res0, "res" = res, "pXY" = difference))
}