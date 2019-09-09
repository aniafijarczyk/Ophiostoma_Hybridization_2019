library(EBImage)
library(dplyr)

Ophiotoma_image2shape_12pos_PixIntens <- function(f){
  
  # Create a folder to store images and results
  
  system(paste0("mkdir ", gsub(".JPG", "", f), "_PixIntens"))
  
  # Read the image
  img = readImage(f)
  #img
  
  # Output the image as read by EBImage with the boundaries that will be croped out
  pdf(paste0(gsub(".JPG", "", f), "_PixIntens/01_CropBoundaries.pdf"))
  display(img, method = "raster")
  abline(v = 584, col = "red")
  abline(v = 4670, col = "red")
  
  abline(h = 450, col = "red")
  abline(h = 3050, col = "red")
  dev.off()
  
  
  # Crop Image to remove the border and inscription on the lid
  img_crop = img[584:4670, 450:3050, 1:3]
  #img_crop
  
  # Output the image croped according to predefined boundaries
  pdf(paste0(gsub(".JPG", "", f), "_PixIntens/02_CropImage.pdf"))
  display(img_crop, method = "raster", all = T)
  dev.off()
  
  # Sum the RGB matrix 
  .list_rgb <- list(imageData(img_crop)[,,1], imageData(img_crop)[,,2], imageData(img_crop)[,,3])
  img_crop_Grayscale <- Reduce("+", .list_rgb)
  
  
  pdf(paste0(gsub(".JPG", "", f), "_PixIntens/03_CropImageGrayscale.pdf"))
  display(img_crop_Grayscale, method = "raster", all = T)
  dev.off()

  dim(img_crop_Grayscale)

  # Boundaries of each column and row
  col1 <- c(125:1025) # 540 pixel wide
  col2 <- c(1025:1925) # 540 pixel wide
  col3 <- c(1925:2825) # 540 pixel wide
  col4 <- c(2825:3725) # 540 pixel wide
  row1 <- c(0:867) # 540 pixel wide
  row2 <- c(867:1734) # 540 pixel wide
  row3 <- c(1734:2601) # 540 pixel wide
  
  
  # Crop each position, output each cropped position and extract the shape information within each position
  
  ##### Position A1
  A1 <- img_crop_Grayscale[col1, row1]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_A1.pdf"))
  #display(A1, method = "raster")
  #dev.off()
  
  A1_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "A1", s.area = sum(A1), version = "pixIntens")
  #####
  
  
  ##### Position A2
  A2 <- img_crop_Grayscale[col2, row1]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_A2.pdf"))
  #display(A2, method = "raster")
  #dev.off()
  
  A2_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "A2", s.area = sum(A2), version = "pixIntens")
  #####
  
  ##### Position A3
  A3 <- img_crop_Grayscale[col3, row1]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_A3.pdf"))
  #display(A3, method = "raster")
  #dev.off()
  
  A3_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "A3", s.area = sum(A3), version = "pixIntens")
  #####
  
  ##### Position A4
  A4 <- img_crop_Grayscale[col4, row1]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_A4.pdf"))
  #display(A4, method = "raster")
  #dev.off()
  
  A4_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "A4", s.area = sum(A4), version = "pixIntens")
  #####
  
  ##### Position B1
  B1 <- img_crop_Grayscale[col1, row2]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_B1.pdf"))
  #display(B1, method = "raster")
  #dev.off()
  
  B1_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "B1", s.area = sum(B1), version = "pixIntens")
  #####
  
  ##### Position B2
  B2 <- img_crop_Grayscale[col2, row2]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_B2.pdf"))
  #display(B2, method = "raster")
  #dev.off()
  
  B2_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "B2", s.area = sum(B2), version = "pixIntens")
  #####
  
  ##### Position B1
  B3 <- img_crop_Grayscale[col3, row2]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_B3.pdf"))
  #display(B3, method = "raster")
  #dev.off()
  
  B3_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "B3", s.area = sum(B3), version = "pixIntens")
  #####
  
  ##### Position B4
  B4 <- img_crop_Grayscale[col4, row2]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_B4.pdf"))
  #display(B4, method = "raster")
  #dev.off()
  
  B4_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "B4", s.area = sum(B4), version = "pixIntens")
  #####
  
  ##### Position C1
  C1 <- img_crop_Grayscale[col1, row3]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_C1.pdf"))
  #display(C1, method = "raster")
  #dev.off()
  
  C1_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "C1", s.area = sum(C1), version = "pixIntens")
  #####
  
  ##### Position C2
  C2 <- img_crop_Grayscale[col2, row3]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_C2.pdf"))
  #display(C2, method = "raster")
  #dev.off()
  
  C2_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "C2", s.area = sum(C2), version = "pixIntens")
  #####
  
  ##### Position C3
  C3 <- img_crop_Grayscale[col3, row3]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_C3.pdf"))
  #display(C3, method = "raster")
  #dev.off()
  
  C3_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "C3", s.area = sum(C3), version = "pixIntens")
  #####
  
  ##### Position C4
  C4 <- img_crop_Grayscale[col4, row3]
  
  #pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_C4.pdf"))
  #display(C4, method = "raster")
  #dev.off()
  
  C4_info <- cbind(picture_path = f, picture_name = gsub(".*/(.*)", "\\1", f), pos = "C4", s.area = sum(C4), version = "pixIntens")
  #####
  
  
  # Grouped all the shape information of each position in a table that is output
  dfs <- list(as.data.frame(A1_info), as.data.frame(A2_info), as.data.frame(A3_info), as.data.frame(A4_info), as.data.frame(B1_info), as.data.frame(B2_info), as.data.frame(B3_info), as.data.frame(B4_info), as.data.frame(C1_info), as.data.frame(C2_info), as.data.frame(C3_info), as.data.frame(C4_info))
  
  info <- Reduce(full_join, dfs)
  
  write.table(info, paste0(gsub(".JPG", "", f), "_PixIntens/14_info.txt"), quote = F, row.names = F)
  
  
}











