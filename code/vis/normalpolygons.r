# Code to make figures with overlapping normal PDFs for schematics
# Uses ggplot2
# Note: the extrafont library will take a long time to load when you first install it. If you skip that, the Helvetica will not display properly.
# After the first time, it will not take any time to load, it just has to import the fonts into R's directory.

# For documentation on stat_function, see http://docs.ggplot2.org/current/stat_function.html
# For a list of all the colors, see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

library(ggplot2)
library(extrafont)

# Set plot parameters
theme_john <- theme_bw() + theme(panel.grid = element_blank(), 
                                 axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 18),
                                 text = element_text(family = 'Helvetica'))

# Add some additional parameters to get rid of the axis ticks and numbers if you want.
theme_noaxisnumbers <- theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

p <- ggplot(data.frame(x = c(-20, 20)), aes(x = x)) + # Creates a dummy plot with x values extending between -20 and 20
  theme_john + # Adds your custom theme
  theme_noaxisnumbers + # Gets rid of axis ticks and numbers, delete if needed
  scale_y_continuous(name = 'The Y axis', limits = c(0, 0.3), expand = c(0, 0)) + # Y axis with name and limits, the expand argument means there is no gap below the 0 value
  scale_x_continuous(name = 'The X axis') +
  stat_function(fun = 'dnorm', geom = 'polygon', fill = 'indianred', alpha = 0.9, args = list(mean = 0, sd = 2)) + # You can make as many of these as you want, hopefully it is obvious how to change the color, mean, and sd. alpha is transparency value where 1 is opaque.
  stat_function(fun = 'dnorm', geom = 'polygon', fill = 'skyblue', alpha = 0.5, args = list(mean = 3, sd = 4)) +
  stat_function(fun = 'dnorm', geom = 'polygon', fill = 'forestgreen', alpha = 0.2, args = list(mean = -4, sd = 2)) + # Change the order of these lines to change the order they appear on the plot.
  geom_segment(x = 10, xend = 10, y = 0, yend = 0.2, color = 'goldenrod', size = 1.5, lineend = 'round') + # Here is a line which you could use for the schematic where there's no variation.
  geom_text(x = 0, y = 0.25, label = 'Hi John', family = 'Helvetica', size = 16) # Example of how to add labels to the plot

# Display plot in the plotting window
p

# Save plot to a png file with high res and 6x6"
ggsave('myplot.png', p, height = 6, width = 6, dpi = 400)

# Save plot to a pdf file
ggsave('myplot.pdf', p, height = 6, width = 6)  
