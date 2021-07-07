source("read_files.R")
source("QTL_plots.R")

exp = 82
plates = c(1)
days = c(4,5,7,10,12,14)
data82 = load_data(exp,plates,days) # returns a list with three dataframes (data, data_end, pam)

exp = 87
plates =c(1,2)
days = c(4,5,7,10,12,14)
data87 = load_data(exp,plates,days)

exp=88
plates=c(1,2,3,4,5)
days = c(4,5,7,10,12,14)
data88 = load_data(exp,plates,days)

exp=92
plates=c(1,2,3,4,5,6)
days = c(4,5,7,10,12,14)
data92 = load_data(exp,plates,days)

exp=93
plates=c(1,2,3,4,5)
days = c(4,5,7,10,12,14)
data93 = load_data(exp,plates,days)

exp=95 
plates=c(1,2,3,4)
days = c(4,5,7,10,12,14)
data95 = load_data(exp,plates,days)

exp=96
plates=c(1,2,3,4,5)
days = c(4,5,7,10,12,14)
data96 = load_data(exp,plates,days)

exp=97
plates=c(1,2,3,4,5)
days = c(4,5,7,10,12,14)
data97 = load_data(exp,plates,days)


exp=98
plates=c(1,2,3,4,5)
days = c(4,5,7,10,12,14)
data98 = load_data(exp,plates,days)

exp=99
plates=c(1,2)
days = c(4,5,7,10,12,14)
data99 = load_data(exp,plates,days)

exp=100
plates=c(1,2)
days = c(4,5,7,10,12,14)
data100 = load_data(exp,plates,days)

exp=103
plates=c(1,2,3,4,5,6,7,8)
days = c(4,5,7,10,12,14)
data103 = load_data(exp,plates,days)

exp=104
plates=c(1,2)
days = c(4,5,7,10,12,14)
data104 = load_data(exp,plates,days)


## Merging experiments data

id = 1
data = merge(data82[[id]],data87[[id]],all = TRUE) %>% merge(data88[[id]], all= TRUE) %>% merge(data92[[id]], all= TRUE) %>% merge(data93[[id]], all= TRUE) %>% merge(data95[[id]], all= TRUE)
data = merge(data,data96[[id]],all = TRUE) %>% merge(data97[[id]], all= TRUE)  %>% merge(data98[[id]], all= TRUE)  %>% merge(data99[[id]], all= TRUE)  %>% merge(data100[[id]], all= TRUE)
data = merge(data,data103[[id]],all = TRUE) %>% merge(data104[[id]], all= TRUE)

id = 2
data_end = merge(data82[[id]],data87[[id]],all = TRUE) %>% merge(data88[[id]], all= TRUE) %>% merge(data92[[id]], all= TRUE) %>% merge(data93[[id]], all= TRUE) %>% merge(data95[[id]], all= TRUE)
data_end = merge(data,data96[[id]],all = TRUE) %>% merge(data97[[id]], all= TRUE) %>% merge(data98[[id]], all= TRUE) %>% merge(data99[[id]], all= TRUE) %>% merge(data100[[id]], all= TRUE)


id = 3
pam1 = merge(data82[[id]],data87[[id]],all = TRUE) %>% merge(data88[[id]], all= TRUE) %>% merge(data92[[id]], all= TRUE) #%>% merge(data93[[id]], all= TRUE, memory.limit(size = 25000))

pam2 = rbind(data82[[id]],data87[[id]])  %>% rbind(data88[[id]]) %>% rbind(data92[[id]])

########################### DATA PLOTS ###########################

strains = unique(data$strain)
strains = strains[strains != "BLANK"]

##Fluorescence plot
get_plots(data=data,strains=strains,colour="green", x_var="day", y_var="chl_fluoro",x_lab="Day",y_lab="RFU",tag="fluoro")

## Nephelometry plot
get_plots(data=data,strains=strains,colour="blue", x_var="day", y_var="neph",x_lab="Day",y_lab="RNU",tag="neph")

## FNc plot
get_plots(data=data,strains=strains,colour="red", x_var="day", y_var="FNc",x_lab="Day",y_lab="RFU/RNU (Control corrected)",tag="FNc")

##Fluoro & Nephelometry plot for BLANKS
get_plots(data=data,strains="BLANK",colour="green", x_var="day", y_var="raw_fluoro",x_lab="Day",y_lab="RFU",tag="fluoro")
get_plots(data=data,strains="BLANK",colour="blue", x_var="day", y_var="raw_neph",x_lab="Day",y_lab="RNU",tag="neph")


### PAM Plots
strains = unique(pam$strain)
strains = strains[strains != "BLANK"]
get_plots(data=pam,strains=strains,colour="orange", x_var="day", y_var="F0",x_lab="Day",y_lab="F0",tag="F0")


# Data Integral calculation

integrand = function(x) 
  predict(lm, newdata = list(day=x)) 

calculate_integral <- function(modelo,col){
  
  modelo$integral={}
  modelo$slope={}
  i = 1
  for (model in modelo$mod){
    lm = model
    value=integrate(integrand,4,14)$value
    modelo$integral[i] = value
    modelo$slope[i] = coef(summary(model))[3]
    i = i + 1
  }
  
  # Adding AUCc
  matrix=modelo %>%  select(strain,replica,plate,exp,integral)
  index=c("exp","plate","strain","replica")
  subcol=subset[subset$day==4, ] %>%  select(strain,replica,plate,exp,col)
  matrix=merge(x=matrix,y=subcol, by.x=index, by.y=index,all = TRUE)
  colnames(matrix) = c(index,"integral","temp")
  matrix$AUCc = matrix$integral/matrix$temp
  
  # Adding model slope
  matrix=merge(x=matrix,y=modelo, by.x=index, by.y=index,all = TRUE)
  matrix=matrix %>% select(exp,plate,strain,replica,integral.y,AUCc,slope)
  auccol=paste("AUCc_",col,sep="")
  slopecol=paste("slope_",col,sep="")
  colnames(matrix) = c(index,"integral",auccol,slopecol)
  return(matrix)
}



#data = data88
subset <- data %>%  filter(type!='B') %>% group_by(strain,replica,plate,exp)
modelo <- subset %>% do(mod = lm(FNc~day + I(day^2), data = .))
matrix <- calculate_integral(modelo,"FNc")
head(matrix,3)


#  Pam integral calculation

subset <- pam %>% group_by(strain,replica,plate,exp)
modelo <- subset %>% do(mod = lm(F0_N~day + I(day^2), data = .))

matrix_tmp <- calculate_integral(modelo,"F0_N")
matrix=merge(x=matrix,y=matrix_tmp, by.x=index, by.y=index,all = TRUE)
matrix=matrix %>% select(index,AUCc_FNc,slope_FNc,WGAc,AUCc_F0_N,slope_F0_N)

head(matrix,3)

#=---------------------------------------------------------------
subset <- pam %>% group_by(strain,replica,plate,exp)
modelo <- subset %>% do(mod = lm(Fm_N~day + I(day^2), data = .))

matrix_tmp <- calculate_integral(modelo,"Fm_N")
matrix=merge(x=matrix,y=matrix_tmp, by.x=index, by.y=index,all = TRUE)
matrix=matrix %>% select(index,AUCc_FNc,slope_FNc,WGAc,AUCc_F0_N,slope_F0_N,AUCc_Fm_N,slope_Fm_N)

head(matrix,3)

#=---------------------------------------------------------------

subset <- pam %>% group_by(strain,replica,plate,exp)
modelo <- subset %>% do(mod = lm(FPQ~day + I(day^2), data = .))

matrix_tmp <- calculate_integral(modelo,"FPQ")
matrix=merge(x=matrix,y=matrix_tmp, by.x=index, by.y=index,all = TRUE)
matrix=matrix %>% select(index,AUCc_FNc,slope_FNc,WGAc,AUCc_F0_N,slope_F0_N,AUCc_Fm_N,slope_Fm_N,AUCc_FPQ,slope_FPQ)

head(matrix,3)

#=---------------------------------------------------------------
subset <- pam %>% group_by(strain,replica,plate,exp)
modelo <- subset %>% do(mod = lm(Fv_Fm~day + I(day^2), data = .))

matrix_tmp <- calculate_integral(modelo,"Fv_Fm")
matrix=merge(x=matrix,y=matrix_tmp, by.x=index, by.y=index,all = TRUE)
matrix=matrix %>% select(index,AUCc_FNc,slope_FNc,WGAc,AUCc_F0_N,slope_F0_N,AUCc_Fm_N,slope_Fm_N,AUCc_FPQ,slope_FPQ,AUCc_Fv_Fm,slope_Fv_Fm)

head(matrix,3)


#=---------------------------------------------------------------


ggplot(matrix, aes(x = strain, y = AUCc)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle=90))

ggplot(matrix, aes(x = strain, y = log (AUCc))) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle=90))

ggplot(matrix_2, aes(x = strain, y = WGAc)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle=90))

ggplot(matrix_2, aes(x = strain, y = log(WGAc))) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle=90))

ggplot(matrix, aes(x = strain, y = slope)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle=90))

subset(matrix, strain == "Ec569-92",)


ggplot(subset(matrix, exp == "93",), aes(x = strain, y = AUCc)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle=90))

ggplot(subset(matrix, exp == "92",), aes(x = strain, y = AUCc)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle=90))

plot(log(matrix$WGAc),log(matrix$AUCc))
plot(matrix$slope,log(matrix$AUCc))

write.csv(matrix,'matrix.csv')
