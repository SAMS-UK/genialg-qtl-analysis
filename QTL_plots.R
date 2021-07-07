


get_plots <- function(data,strains,colour,x_var, y_var,x_lab,y_lab,tag){
  for (strain in strains){
    datos = data[data$strain==strain,]
    inter <- "interaction(exp,plate,replica)"
    plot = ggplot(datos, aes_string(x = x_var, y = y_var,
                             group = inter,
                             colour="strain")) + 
      geom_point(colour = colour)+
      geom_line(colour = colour) +
      ylim(c(0,max(na.omit(datos %>% select(y_var))))) +
      xlab(x_lab) + ylab(y_lab)+
      theme_bw() +
      guides(fill=FALSE)
    print(strain)
    image_path = paste("images/",tag,"_",strain,".png",sep="")
    print(image_path)
    ggsave(image_path,plot,width = 5, height = 4, dpi = 300, device='png')
  }
}
