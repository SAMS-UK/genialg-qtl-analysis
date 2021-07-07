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


## Merging experiments data

id = 1
data = merge(data82[[id]],data87[[id]],all = TRUE) %>% merge(data88[[id]], all= TRUE) %>% merge(data92[[id]], all= TRUE) %>% merge(data93[[id]], all= TRUE) %>% merge(data95[[id]], all= TRUE)

id = 2
data_end = merge(data82[[id]],data87[[id]],all = TRUE) %>% merge(data88[[id]], all= TRUE) %>% merge(data92[[id]], all= TRUE) %>% merge(data93[[id]], all= TRUE) %>% merge(data95[[id]], all= TRUE)

id = 3
pam = merge(data82[[id]],data87[[id]],all = TRUE) %>% merge(data88[[id]], all= TRUE) %>% merge(data92[[id]], all= TRUE) #%>% merge(data93[[id]], all= TRUE, memory.limit(size = 25000))


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


