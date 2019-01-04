a = 1
sigmaSq = 1
confLevel = 0.99 #gamma
alpha = 0.05
alphaZ = 2
sigmaSqZ = 1
pZ = 0.6
k = 0.3

N = 100

source("basicMethods.R")
source("Quantiles.R")

smpl1 = rnorm(N, a, sigmaSq^(1/2))
smpl2 = rnorm(N, a, sigmaSq^(1/2))

#Таблицы выборок
smpl1Str = list("Выборка" = paste(smpl1, collapse = ", "))
smpl1Table = flextable(as.data.frame(smpl1Str))
smpl1Table = width(smpl1Table, width = 5.7)
smpl2Str = list("Выборка" = paste(smpl2, collapse = ", "))
smpl2Table = flextable(as.data.frame(smpl2Str))
smpl2Table = width(smpl2Table, width = 5.7)
#------------------------------------------------------------

paramsTable = flextable(as.data.frame(list("a" = a, "sigmaSq" = sigmaSq, "gamma" = confLevel, "alpha" = alpha,
                  "alpha0" = alphaZ, "sigmaSq0" = sigmaSqZ, 
                  "p0" = pZ, "k" = k)))

#Описание выборок
desc = as.data.frame(getSampleDescription(smpl1))
descNames = make.names(c("Выборочное среднее", "Исправленная выборочная дисперсия", "Выборочная дисперсия",
              "Медиана", "Минимум", "Максимум", "Размах", "Коэффициент вариации", "Коэффициент осциляции"))
colnames(desc) = descNames
tableSmpl1 = flextable(desc, col_keys = descNames)
desc = as.data.frame(getSampleDescription(smpl2))
colnames(desc) = descNames
tableSmpl2 = flextable(desc, col_keys = descNames)
#--------------------------------------------------------------

#Гистограмма и функция распределения
par(mar = rep(2, 4))
histogramm = function() print(hist(smpl1, breaks = 6, freq = FALSE))
png(filename = "Гистограмма.png")
hist(smpl1, breaks = 6, freq = FALSE, main = "Гистограмма с плотностью", xlab = "")
lines(density(smpl1), col = "red")
dev.off()
png(filename = "Распределение.png")
plot(ecdf(smpl1), main = "Эмпирическая функция распределения")
dev.off()
#--------------------------------------------------------------

#Доверительный интервал при известной дисперсии
confidenceIntervalEV = confidenceIntervalEV(smpl1, sigmaSq, confLevel)
#---------------------------------------------

#Доверительный интервал при неизвестной дисперсии
confidenceIntervalEVFull = confidenceIntervalEVFull(smpl1, confLevel)
#-----------------------------------------------------------------

#Доверительный интервал для дисперсии
confidenceIntervalDisp = confidenceIntervalDisp(smpl1, confLevel)
#------------------------------------

#Согласованность с нормальным распределением
ndCheck = normalDistributionCheck2(smpl1, alpha)
ndTable = flextable(ndCheck$df)
#------------------------------------------------------------------

#Однородность
hmCheck = homogeneityCheck(smpl1, smpl2, alpha)
hmTable = flextable(hmCheck$df)
#------------------------------------------------

#Гипотезы о параметрах
paramsH = dispAndEVCheck(smpl1, alphaZ, sigmaSqZ, alpha)
meansEq = meansEqualCheck(smpl1, smpl2, alpha)
dispEq = dispEqualCheck(smpl1, smpl2, alpha)
#------------------------------------------------

#Анализ корреляции
corrCheck = corellationCheck(smpl1, smpl2, confLevel, alpha)
#--------------------------------------------------

#Анализ вероятностей
probCheck = probAnalysis(smpl1, smpl2, pZ, a, k, (sigmaSq)^(1/2), confidenceLevel = confLevel, alpha = alpha)
#-------------------------------------------------------


tables = list(params = paramsTable, tableSmpl1 = tableSmpl1,
              tableSmpl2 = tableSmpl2, smpl1 = smpl1Table,
              smpl2 = smpl2Table, 
              ndCheck = ndTable, hmCheck = hmTable)
plots = list(histogramm = histogramm)

print("Generating tables")
body_add_flextables("CourseTemplate.docx","CourseTemplateOut.docx", tables)
print("Tables generated. Generating inline R code")
renderInlineCode("CourseTemplateOut.docx","CourseTemplateOut.docx")
print("R code generated")
print("Don't forget to check document, insert your name, pictures and create different design. Also check for NULL values and errors. Re-gen doc if needed")