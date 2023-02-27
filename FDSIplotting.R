library(ggplot2)
library(raster)
library(raster)
library(ncdf4)
library(ggplot2)
library(magick)

color_pal <- colorRampPalette(c("blue", "lightblue",
                                "white", "orange","red"))


fdsi36=brick("C:/Users/Brown/Documents/FLASH/FDSI.nc")
fdsi9=brick("C:/Users/Brown/Documents/FLASH/9kmV2/FDSI30_3.nc")
tsOut=as.numeric(extract(fdsi9, cbind(-115,40)))
plot(tsOut, type="l")

fdsi36crop = crop(fdsi36, fdsi9)

plot(fdsi36crop[[4]])
plot(fdsi9[[4]])

brk=c(0,0.2,0.4,0.5,0.6,0.71,1)


#convert rasterbrick to dataframe
fdsi36data = as.data.frame(as(fdsi36crop, "SpatialPixelsDataFrame")) # Spatial dataframe of SM data)

fdsi9data = as.data.frame(as(fdsi9, "SpatialPixelsDataFrame")) # Spatial dataframe of SM data)

value = fdsi9data[[1030]]
value = fdsi36data[[1030]]
{
  globalPlot=ggplot() +                             # Initialize the plot
    geom_tile(data=fdsi36data,                    # Add soil moisture data
              aes_string(x="x", y="y", fill= cut(value, breaks=brk)),alpha=1)+       # Provide X, Y and Z data
    scale_fill_manual(drop=TRUE, values=(color_pal(length(brk)-1)), 
                      na.value="transparent", name="FDSI", 
                      labels=levels(cut(value, breaks=brk)))+           # Use user-defined colormap                         # Set limits of colorbar
                    coord_cartesian(xlim = c(-104.9998, -89.9624),
                    ylim = c(25.0345, 50.0269))  +                                            #xlab("Longitude")+                            # X-axis label
    ylab("Latitude")+  xlab("Longitude")+                            # Y-axis label
    theme_bw() +
    borders("world",                              # Add global landmass boundaries
            colour="gray43",  size = 0.1,                    # Fill light-gray color to the landmass
            fill="transparent")  +                # Transparent background 
    borders("state",                # Add US state borders
            colour = "gray43",      # Use light-gray color
            fill = "transparent")+
    labs(title = 'Central US Drought Outlook - FLASH 1.0 - 36km',
         subtitle = "January 24, 2018",
         caption = 'FDSI >0.71 indicates emerging/ sustained flash droughts')
  
  print(globalPlot)
}  
ggsave(plot = globalPlot, path = "C:/Users/brown/Documents/", filename = paste("flash1v2_36km_R3_01-24-2018.png"), bg="white", units="in", dpi=300, width=5, height=5)



hist(fdsi9data[[1030]], xlim=c(0,1), xlab="FDSI Values", main = "FLASH 2.0 - 9km, January 24, 2018")

hist(fdsi36data[[1030]], xlim=c(0,1), xlab="FDSI Values", main = "FLASH 1.0 - 36km, January 24, 2018")
