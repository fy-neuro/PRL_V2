library(pdftools)
library(magick)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/fanyi/code/")

merge_pdfs_horizontally <- function(pdf1, pdf2, output_pdf) {
  # Read the pages of the two PDFs as images
  pdf1_images <- image_read_pdf(pdf1)
  pdf2_images <- image_read_pdf(pdf2)
  
  # Check if the two PDFs have the same number of pages
  if (length(pdf1_images) != length(pdf2_images)) {
    stop("The two PDFs have different numbers of pages. Cannot merge horizontally!")
  }
  
  # Create a list to store horizontally merged pages
  combined_images <- vector("list", length(pdf1_images))
  
  # Loop through each page and combine them horizontally
  for (i in seq_along(pdf1_images)) {
    combined_images[[i]] <- image_append(c(pdf1_images[i], pdf2_images[i]), stack = FALSE)
  }
  
  # Combine all pages into a single PDF
  combined_pdf <- image_write_pdf(image = image_join(combined_images), path = output_pdf)
  
  # Return the output path for confirmation
  return(combined_pdf)
}

# Example Usage
merge_pdfs_horizontally("../plots/all_pathways_plots_Control_1.pdf", "../plots/all_pathways_plots_PD_1.pdf", "output_combined.pdf")

