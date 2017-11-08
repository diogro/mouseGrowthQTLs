library(geomorph)
library(sp)

raw_p_data = read_csv("./data/Mandibles/parental_strains_landmarks.csv")

landmarks = geomorph::arrayspecs(A = raw_p_data[,5:34], p = 15, k = 2, sep = "")
dimnames(landmarks)[[3]] = raw_p_data$ID

plot(-landmarks[,,2])
text(-landmarks[,,2], labels=1:15, cex= 3.7)

area_borders = list(
  area1 = c(4, 5, 6),
  area2 = c(6, 7, 8, 9),
  area3 = c(9, 10, 11, 12),
  area4 = c(4, 6, 9, 12),
  area5 = c(2, 4, 12, 13),
  area6 = c(2, 3, 4),
  area7 = c(1, 2, 13, 14, 15)
)

get_region_polygon = function(area_lands, ind){
  box.coords <- ind[area_lands,]
  box.hpts <- chull(x = box.coords[,1], y = box.coords[,2])
  box.hpts <- c(box.hpts, box.hpts[1])
  box.chull.coords <- box.coords[box.hpts,]
  chull.poly <- Polygon(box.chull.coords, hole=F)
  chull.poly
}
get_region_areas = function(ind){
  llply(area_borders, get_region_polygon, ind) %>% laply(function(x) x@area)
}
parental_areas = t(apply(landmarks, 3, get_region_areas))
colnames(parental_areas) <- area_traits
ancestral_areas = ddply(data.frame(strain = raw_p_data$Strain, parental_areas), .(strain), numcolwise(mean))
# 
# 
# llply(area_borders, get_region_polygon, landmarks[,,1])
# 
# Plot_ConvexHull<-function(area, lcolor){
#   xcoord = area[,1]
#   ycoord = area[,2]
#   hpts <- chull(x = area[,1], y = area[,2])
#   hpts <- c(hpts, hpts[1])
#   lines(xcoord[hpts], ycoord[hpts], col = lcolor)
# }  
# ind = landmarks[,,2]
# xrange <- range(ind[,1])
# yrange <- range(ind[,2])
# plot(ind[area_borders[[i]],], type = "p", pch = 1, col = "black", xlim = c(xrange), ylim = c(yrange))
# Plot_ConvexHull(ind[area_borders[[i]],], lcolor = "black")
# for(i in seq_along(area_borders)){
#   points(ind[area_borders[[i]],], type = "p", pch = 1, col = "black", xlim = c(xrange), ylim = c(yrange))
#   Plot_ConvexHull(ind[area_borders[[i]],], lcolor = "black")
# }
