require(ggplot2, quietly = TRUE)

# Modify ggplot2 defaults -------------------------------------------------

# This script will check for the following variables in the current environment:
# base.size (ggplot2 font size)
# fig.dir   (directory in which figures will be saved)
# prefix    (text appended to files generated with save_plot)

if(!exists("base_size")) {
  base_size <- 9
}

if(!exists("fig.dir")) {
  fig.dir <- getwd()
}

if(!exists("fig.prefix")) {
  fig.prefix <- ""
}


# Update grey theme to remove margins and use set font size to base_size
theme_set(theme_grey(base_size))
theme_update(plot.margin = rep(grid::unit(0, "lines"), 4))


# Convenience function for saving figures ---------------------------------

save_plot <- function(filename, width = 5, height = 4, path = fig.dir, 
  dpi = 150, prefix = fig.prefix, ...) {
  
  filename <- paste(prefix, filename, sep = "")
  
  ggsave(filename, width = width, height = height, path = path, 
         dpi = dpi, ...)
}
