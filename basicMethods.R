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

  hist(testX, breaks = steps, freq = FALSE)
  lines(density(testX)) #плотность
  plot(ecdf(testX)) #эмпирическая функция распределения
}