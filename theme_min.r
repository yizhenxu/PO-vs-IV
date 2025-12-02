theme_min <- function(size = 12, font = "sans", textColor = 'black', whiteBG = TRUE, axisSize = 0.5,
											labsize = NA, ...) {

	if (!whiteBG) {
		if (textColor=='black') {
			textColor='white'
		}
	}

	if (is.na(labsize)) {
		labsize <- size
	}

  element_text = function(...)
    ggplot2::element_text(family=font, colour=textColor, ...)
	#element_line = function(...)
	#	ggplot2::element_line(colour=textColor,...)

	theme(
    axis.title.x = element_text(face="bold",size=size),
    axis.title.y = element_text(face="bold",angle=90,size=size),
    axis.text.x = element_text(size=labsize),
    axis.text.y = element_text(size=labsize),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), # no box, axis labels instead
    legend.background = element_blank(),
    legend.key = element_blank(), # no box around legend
    legend.title = element_blank(), # no legend title
    legend.text = element_text(size=size),
    panel.background = element_blank(),
    panel.grid = element_line(colour="grey92"),
    panel.grid.minor = element_line(linewidth=rel(0.5)),
		plot.background = element_blank(),
    strip.background=element_blank(), # no box on panel labels
    strip.text.x = element_text(hjust=0,size=size),
    strip.text.y = element_text(size=size,angle=0) # make row panel labels horizontal
		#,size=size) # dropping this with the switch from opts to theme probably broke something
	)
}
