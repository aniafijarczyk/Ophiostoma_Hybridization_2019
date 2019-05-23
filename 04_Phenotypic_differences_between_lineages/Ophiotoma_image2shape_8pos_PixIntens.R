#!/usr/bin/env Rscript


#
# Image analyses for Ophiostoma 
#
# usage:
# Ophiotoma_image2shape_8pos f=~/Dropbox/Postdoc/Ophiostoma/Test_2_Robot_batch73/d000073_220_002_18-05-24_11-46-20.JPG
#
#



library(EBImage)
library(dplyr)

#argz <- commandArgs()

# f = "/Users/Lenou/Dropbox/Postdoc/Ophiostoma/Test_2_Robot_batch73/d000073_220_001_18-05-23_11-45-37.JPG"

#f = "~/Dropbox/Postdoc/Ophiostoma/Test_2_Robot_batch73/d000073_220_002_18-05-24_11-46-20.JPG"

Ophiotoma_image2shape_8pos_PixIntens <- function(f){

# Create a folder to store images and results
system(paste0("mkdir ", gsub(".JPG", "", f), "_PixIntens"))

# Read the image
img = readImage(f)
#img

# Output the image as read by EBImage with the boundaries that will be croped out
pdf(paste0(gsub(".JPG", "", f), "_PixIntens/01_CropBoundaries.pdf"))
display(img, method = "raster")
abline(v = 484, col = "red")
abline(v = 4700, col = "red")
abline(h = 400, col = "red")
abline(h = 3056, col = "red")
dev.off()


# Crop Image to remove the border and inscription on the lid
img_crop = img[484:4700, 400:3056, 1:3]
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

rect(10, 0, 550, 540, col = "red", density = 0)
rect(10, 2117, 550, 2657, col = "red", density = 0)
rect(1360, 0, 1900, 540, col = "red", density = 0)
rect(1360, 2117, 1900, 2657, col = "red", density = 0)
rect(2417, 0, 2957, 540, col = "red", density = 0)
rect(2417, 2117, 2957, 2657, col = "red", density = 0)
rect(3677, 0, 4217, 540, col = "red", density = 0)
rect(3677, 2117, 4217, 2657, col = "red", density = 0)

dev.off()




# Boundaries of each column and row
col1 <- c(10:550) # 540 pixel wide
col2 <- c(1360:1900) # 540 pixel wide
col3 <- c(2417:2957) # 540 pixel wide
col4 <- c(3677:4217) # 540 pixel wide
row1 <- c(0:540) # 540 pixel wide
row2 <- c(2117:2657) # 540 pixel wide


# Crop each position, output each cropped position and extract the shape information within each position

##### Position A1
A1 <- img_crop_Grayscale[col1, row1]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/06_A1.pdf"))
display(A1, method = "raster")
dev.off()

A1_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "A1", s.area = sum(A1), version = "pixIntens")
#####


##### Position A8
A8 <- img_crop_Grayscale[col1, row2]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/07_A8.pdf"))
display(A8, method = "raster")
dev.off()

A8_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "A8", s.area = sum(A8), version = "pixIntens")
#####


##### Position B1
B1 <- img_crop_Grayscale[col2, row1]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/08_B1.pdf"))
display(B1, method = "raster")
dev.off()

B1_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "B1", s.area = sum(B1), version = "pixIntens")
#####


##### Position B8
B8 <- img_crop_Grayscale[col2, row2]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/09_B8.pdf"))
display(B8, method = "raster")
dev.off()

B8_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "B8", s.area = sum(B8), version = "pixIntens")
#####


##### Position C1
C1 <- img_crop_Grayscale[col3, row1]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/10_C1.pdf"))
display(C1, method = "raster")
dev.off()

C1_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "C1", s.area = sum(C1), version = "pixIntens")
#####


##### Position C8
C8 <- img_crop_Grayscale[col3, row2]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/11_C8.pdf"))
display(C8, method = "raster")
dev.off()

C8_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "C8", s.area = sum(C8), version = "pixIntens")
#####


##### Position D1
D1 <- img_crop_Grayscale[col4, row1]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/12_D1.pdf"))
display(D1, method = "raster")
dev.off()

D1_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "D1", s.area = sum(D1), version = "pixIntens")
#####


##### Position D8
D8 <- img_crop_Grayscale[col4, row2]

pdf(paste0(gsub(".JPG", "", f), "_PixIntens/13_D8.pdf"))
display(D8, method = "raster")
dev.off()

D8_info <- cbind(picutre_path = f, picutre_name = gsub(".*/(.*)", "\\1", f), pos = "D8", s.area = sum(D8), version = "pixIntens")
#####






# Grouped all the shape information of each position in a table that is output
dfs <- list(as.data.frame(A1_info), as.data.frame(A8_info), as.data.frame(B1_info), as.data.frame(B8_info), as.data.frame(C1_info), as.data.frame(C8_info), as.data.frame(D1_info), as.data.frame(D8_info))

info <- Reduce(full_join, dfs)

write.table(info, paste0(gsub(".JPG", "", f), "_PixIntens/14_info.txt"), quote = F, row.names = F)


}











