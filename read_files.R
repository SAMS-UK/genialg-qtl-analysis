library (ggplot2)
library(readxl)
#library(compare)
library(repr)
library(AUC)
library(stringr)
library(roxygen2)
library(dplyr)


load_data <- function(exp,plates,days){
    raw_data <- read_experiments(exp,plates,days)
    data <- data_processing(raw_data$data)
    data_end = data_processing(raw_data$endpoint)
    pam <- read_experiments_pam(exp,plates,days)
    pam <- pam_data_processing(pam,data)
    return(list(data,data_end,pam))
}



#' Read fluorescence and nephelometry data in excel files.
#' 
#' @param exp Experiment number.
#' @param plates vector of plates used in the experiment.
#' @param days vector of days used in the experiment.
#' @return dataframe that combines fluoresence and nephelometry data

read_experiments <- function(exp,plates,days){
    
    first = TRUE
    firstendpoint = TRUE
    for(plate in plates){
        for (day in days){
            # Read Fluoro excel file and skip first 12 rows
            filename = paste('data/exp0',exp,'/fluoro_',plate,'_',day,'d.xlsx',sep="")
            #print(filename)
            fluoro=read_excel(filename,skip=12)
            # Filtering useful columns
            fluoro = fluoro[4:11]
            # Renaming columns
            names(fluoro)<- c("raw_fluoro","chl_fluoro", "exp","plate","strain","type","replica","day")

            
            # Read Nephelo excel file
            filename = paste('data/exp0',exp,'/nephelo_',plate,'_',day,'d.xlsx',sep="")
            nephelo=read_excel(filename,skip=12)
            nephelo=nephelo[4:11]
            names(nephelo)<- c("raw_neph","neph", "exp","plate","strain","type","replica","day")
            
            # Merging Fluoro and Nephelo dataframes for a day and a plate
            index = c("exp","plate","strain","type","replica","day")
            merged = merge(x=fluoro,y=nephelo, by.x=index, by.y=index,all = TRUE)            
            
            # Merging merged datasets for all days and plates
            if (first){
                data = merged
                first =FALSE
            }else{
                data = merge(data,merged,all = TRUE)
            }
        }
        # Reading Fluoro Endpoints
        filename = paste('data/exp0',exp,'/fluoro_',plate,'_14d_endpoint.xlsx',sep="")
        #print(filename)
        fluoro=read_excel(filename,skip=12)
        
        # Filtering useful columns
        fluoro=fluoro %>% select(4,6,7,9,10:15)
        
        # Renaming columns
        names(fluoro)<- c("raw_chl_fluoro","raw_fluoro_WGA","chl_fluoro","fluoro_WGA", "exp","plate","strain","type","replica","day")
        
        # Reading Neph Endpoints
        filename = paste('data/exp0',exp,'/nephelo_',plate,'_14d_endpoint.xlsx',sep="")
        
        nephelo=read_excel(filename,skip=12)
        
        # Filtering useful columns
        nephelo=nephelo %>% select(4:11)
        
        # Renaming columns
        names(nephelo)<- c("raw_neph","neph", "exp","plate","strain","type","replica","day")
        
        # Merging Fluoro and Nephelo dataframes for a day and a plate
        index = c("exp","plate","strain","type","replica","day")
        endpoint = merge(x=fluoro,y=nephelo, by.x=index, by.y=index,all = TRUE)
        if (firstendpoint){
            data_end = endpoint
            firstendpoint =FALSE
        }else{
            data_end = merge(data_end,endpoint,all = TRUE)
        }
 
        
    }
    strain_err = data %>% filter(is.na(strain))
    if (nrow(strain_err) > 0) {print(strain_err)}
    data <- list("data" = data, "endpoint" = data_end)
    return(data)
}



read_experiments_pam <- function(exp,plates,days){
    first = TRUE
    for (plate in plates){
        for (day in days){
            filename = paste('data/exp0',exp,'/0',exp,'_',plate,"_",day,'D.RES',sep="")
            print(filename)
            #classes=c("character","numeric","numeric","numeric","numeric","numeric","numeric","character","character","numeric","numeric")
            pam=read.table(filename,skip=0,sep = "\t",header = TRUE,dec = ".",stringsAsFactors = FALSE)[-1,]
            print(nrow(pam))
            pam = pam[c("Fo","Fm","Fv.Fm","Sigma","Experiment","Plate","Strain","Type","Replica","Day")]
            names(pam)<- c("F0","Fm","Fv_Fm", "Sigma","exp","plate","strain","type","replica","day")
            #pam = pam %>% filter(F0!="---------") # Second line in the file contains that symbols
            pam = transform(pam, F0 = as.numeric(F0),
                            Fm = as.numeric(Fm),
                            Fv_Fm = as.numeric(Fv_Fm),
                            Sigma = as.numeric(Sigma))
            if (first){
                    data = pam
                    first =FALSE
                }else{
                    #data = merge(data,pam,all = TRUE) # Inefficient merge
                    data = rbind(data,pam)
            }
            }
    }
    return(data)
}


#' Calculates Control value, ratio_F, ratio_N and FNc from the fluoro and neph merged data
#'
#' @param data dataframe containing fluoro and neph values.
#' @return dataframe with calculated columns Ctrl_f,Ctrl_n,ratio_F, ratio_N and FNc


data_processing <- function(data){
    # Creating new two columns to link fluoro and neph controls 'C' to 'A' (Anisolpidium-challenged) type
    data$Ctrl_f = data$chl_fluoro
    data$Ctrl_n = data$neph

    # Calculating total number of rows
    nrows = nrow(data)

    # Assigning mean Control value to their respective Anisolpidium treatment
    for (i in 1:nrows){
        row = data[i,]
        
        if (row$type == 'A'){
            # Finding Control values for matching plate, exp, strain and day
            query = data[data$exp==row$exp & data$plate==row$plate & data$strain==row$strain & data$day==row$day & data$type=='C',]
            
            # Checking Endpoint case
            if("fluoro_WGA" %in% colnames(data)){
                f = query$fluoro_WGA
            }else{
                f = query$chl_fluoro
            }
            
            Ctrl_f = mean(f)
            Ctrl_n = mean(query$neph)
            
            # Updating Control value in the data dataframe
            data[i,]$Ctrl_f = Ctrl_f
            data[i,]$Ctrl_n = Ctrl_n
        }
    }

    # Removing Control data
    data = data[data$type!='C',]
    
    # Calculating columns ratio_F, ratio_N and FNc
    # Checking Endpoint case
    data$ratio_N = data$neph/data$Ctrl_n
    if("fluoro_WGA" %in% colnames(data)){
        data$ratio_WGA = data$fluoro_WGA/data$Ctrl_f
        data$WGAc = data$ratio_WGA/data$ratio_N
        data %>% rename(Ctrl_WGA = Ctrl_f)
    }else{
        data$ratio_F = data$chl_fluoro/data$Ctrl_f
        data$FNc = data$ratio_F/data$ratio_N
    }
    
    return(data)
}


#' Calculates Control value, ratio_F0, ratio_Fm, F0/ratio_N and Fm/ratio_N
#'
#' @param pam 
#' @param data 
#' @return pam dataframe with calculated columns ratio_F0, ratio_Fm, F0/ratio_N and Fm/ratio_N


pam_data_processing <- function(pam,data){
    data <- data %>%  filter(type!='B')
    # Creating new two columns 
    pam$Ctrl_F0 = pam$F0
    pam$Ctrl_Fm = pam$Fm

    # Calculating total number of rows
    nrows = nrow(pam)

    # Assigning mean Control value to their respective Anisolpidium treatment
    for (i in 1:nrows){
        row = pam[i,]
        
        if (row$type == 'A'){
            # Finding Control values for matching plate, exp, strain and day
            query = pam[pam$exp==row$exp & pam$plate==row$plate & pam$strain==row$strain & pam$day==row$day & pam$type=='C',]
            Ctrl_F0 = mean(query$F0)
            Ctrl_Fm = mean(query$Fm)
            
            # Updating Control value in the data dataframe
            pam[i,]$Ctrl_F0 = Ctrl_F0
            pam[i,]$Ctrl_Fm = Ctrl_Fm
        }
    }

    # Removing Control data
    pam = pam[pam$type!='C',]
    
    #Calculating ratio_F0, ratio_Fm
    pam$ratio_F0 = pam$F0/pam$Ctrl_F0
    pam$ratio_Fm = pam$Fm/pam$Ctrl_Fm
    
    # Calculating F0/ratio_N and Fm/ratio_N
    index = c("exp","plate","strain","type","replica","day")
    pam = merge(x=pam,y=data, by.x=index, by.y=index,all = TRUE)
    
    pam$F0_N=pam$ratio_F0/pam$ratio_N
    pam$Fm_N=pam$ratio_Fm/pam$ratio_N
    pam<-pam[c(index,"Fv_Fm","Sigma","F0","Ctrl_F0","ratio_F0","Fm","Ctrl_Fm","ratio_Fm","ratio_N","F0_N","Fm_N")]
    return(pam)
}

