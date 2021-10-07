# Test script for running R from terminal
install.packages("ggplot2")

min.dataframe <- data.frame(a = 1:10, b = 1:10*3) 

folder_name = 'Output_data'

save(min.dataframe, file=paste(folder_name,"/test_dataframe.Rda", sep=""))

test_plot <- ggplot(data = min.dataframe) + geom_point(aes(x = a, y = b))

ggsave('test_plot.png',
       plot = test_plot,
       device = "png",
       path = folder_name,
       height = 5, width = 8, 
       dpi = "retina"
)

ggsave('test_plot.pdf',
       plot = test_plot,
       device = "pdf",
       path = folder_name,
       height = 5, width = 8, 
       dpi = "retina"
)
