my_theme = theme(
  panel.background = element_blank(),
  plot.background = element_blank(), 
  legend.background = element_rect(fill="transparent", colour=NA),
  legend.key = element_rect(fill="transparent", colour=NA),
  plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
  title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
  axis.text.y = element_text(size = 10, margin = margin(r = 5),hjust=1,vjust=0.5, family="Arial", face="bold",colour="black"),
  axis.text.x = element_text(size = 10, margin = margin(t = 5),hjust=0.5,vjust=1, family="Arial", face="bold",colour="black"), 
  axis.title.y = element_text(size = 12, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, family="Arial", face="bold"),
  axis.title.x = element_text(size = 12, margin = margin(t = 10),hjust=0.5,vjust=1, family="Arial", face="bold"),
  legend.text=element_text(size=10, family="Arial"),
  legend.title=element_text(size=12,family ="Arial"), 
  legend.key.size=unit(1,"line")
)
