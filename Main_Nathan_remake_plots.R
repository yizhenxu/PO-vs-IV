
##########################
## Cumulative risk PDFs ##
##########################

# ---- Packages ----
library(ggplot2)
library(data.table)
library(patchwork)

# ---- Prep factors & colors ----
# facet labels
hazdt[, faclab := factor(analysis,
                         levels = c("ITT","per-protocol","as-treated"),
                         labels  = c("Intention-to-Treat","Receipt of Treatment - Compliers","Receipt of Treatment"))]
# linetype (FALSE = Unadjusted, TRUE = Adjusted)
hazdt[, weighted_f := factor(weighted, levels = c(FALSE, TRUE),
                             labels = c("Unadjusted","Adjusted"))]

# colors for IV/PO (fallback if you don't already have plotColors)
if (!exists("plotColors")) {
  plotColors <- c("IV"="#5FA55A", "PO"="#5672D8")
}

# ---- LEFT: IV vs PO densities (no legend) ----
pLeft <- ggplot(
  hazdt[trtmnt %in% c("PO","IV") & cutpoint == 365 &
          analysis %in% c("ITT","per-protocol","as-treated")]
) +
  geom_density(
    aes(x = cumrisk, colour = trtmnt, linetype = weighted_f),
    linewidth = 0.9, key_glyph = "path"     # <-- legend keys draw as lines
  ) +
  facet_wrap(vars(faclab), ncol = 1, strip.position = "right") +
  scale_color_manual(values = plotColors) +
  scale_linetype_manual(values = c("Unadjusted" = "solid", "Adjusted" = "22")) +
  scale_x_continuous("1-year risk of recurrence", labels = scales::percent) +
  ylab("Density") +
  theme_min(base_family = "sans") +
  theme(legend.position = "none") +
  # optional in-plot labels for IV/PO (remove if not desired)
  geom_text(
    data = data.table(
      x = c(0.23, 0.23, 0.23, 0.50, 0.50, 0.50),
      y = rep(c(8, 8, 8), 2),
      label = rep(c("IV","PO"), each = 3)
    ),
    aes(x = x, y = y, label = label, colour = label),
    hjust = 1, family = "sans"
  )+
  guides(linetype = "none", colour = "none", fill = "none") +
  theme(legend.position = "none") 

# ---- RIGHT: PO - IV densities (provide the legend here) ----
subDT <- hazdt[trtmnt == "PO-IV" & cutpoint == 365 &
                 analysis %in% c("ITT","per-protocol","as-treated")]

library(grid)  # for unit()

pRight <- ggplot(rightDT,
                 aes(x = time, y = value, linetype = adj_f)) +
  # confidence ribbons per adj_f (no legend for alpha/fill)
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = adj_f),
              alpha = 0.18, colour = NA, show.legend = FALSE) +
  # lines on top
  geom_line(colour = "black", linewidth = 1, key_glyph = "path") +
  facet_wrap(vars(faclab), ncol = 1) +
  scale_linetype_manual(values = c("Unadjusted" = "solid", "Adjusted" = "22"),
                        name = NULL) +
  scale_fill_manual(values = c("Unadjusted" = "grey40", "Adjusted" = "grey20"),
                    guide = "none") +
  guides(
    linetype = guide_legend(
      override.aes = list(colour = "black", linewidth = 1.2),
      keywidth = unit(22, "pt"), keyheight = unit(6, "pt")
    )
  ) +
  scale_x_continuous(expression(bold("Days"))) +
  scale_y_continuous(position = "right") +
  theme_min(base_family = "sans") +
  theme(axis.title.y = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center")+
  annotate("text", label = "PO \u2212 IV", x = Inf, y = Inf,
           hjust = 1.05, vjust = 1.3)


# ---- Combine panels and collect the single legend (from pRight) ----
final <- (pLeft |
            pRight) +
  plot_layout(widths = c(7, 5.4), guides = "collect") &
  theme(legend.position = "bottom",
        legend.justification = "center")

# Draw
print(final)

if (writeOutput) {
  ggplot2::ggsave(
    filename = "cumriskPDFs.pdf",
    plot     = final,
    width    = 12.4,   # inches
    height   = 7,      # inches
    units    = "in",
    device   = grDevices::cairo_pdf
  )
  #ggsave(final, file="cumriskPDFs.svg")
  #quartz.save(file='cumriskPDFs.pdf', type="pdf")
}

##############################
## Cumulative Distrn curves ##
##############################
cdfdt[, `:=`(adj=factor(fifelse(weighted, "adjusted", "unadjusted"), levels=c("unadjusted", "adjusted")), analysis=factor(analysis, levels=c("ITT", "as-treated", "per-protocol")))]
cdfdt[, faclab := factor(analysis,
                         levels = c("ITT","per-protocol","as-treated"),
                         labels  = c("Intention-to-Treat","Receipt of Treatment - Compliers","Receipt of Treatment"))]
# linetype (FALSE = Unadjusted, TRUE = Adjusted)
cdfdt[, weighted_f := factor(weighted, levels = c(FALSE, TRUE),
                             labels = c("Unadjusted","Adjusted"))]


library(ggplot2)
library(data.table)
library(patchwork)
library(scales)

# ensure data.table
plotdt <- as.data.table(cdfdt)

# Make a clean "Adjusted/Unadjusted" factor if needed
if (!"adj_f" %in% names(plotdt)) {
  plotdt[, adj_f := ifelse(isTRUE(weighted), "Adjusted", "Unadjusted")]
}
plotdt[, adj_f := factor(adj_f, levels = c("Unadjusted","Adjusted"))]
plotdt$group[plotdt$group=="Oral"] = "PO"

# Optional: colors for treatment groups
plotColors <- c("IV"="#5FA55A", "PO"="#5672D8")

# -------- LEFT: IV vs PO (no ribbons) --------
leftDT <- plotdt[group %in% c("IV","PO") &
                   analysis %in% c("ITT","per-protocol","as-treated")]

pLeft <- ggplot(leftDT,
                aes(x = time, y = value,
                    colour = group, linetype = weighted_f)) +
  geom_line(linewidth = 0.9, key_glyph = "path") +
  facet_wrap(vars(faclab), ncol = 1, strip.position = "right") +
  scale_color_manual(values = plotColors, name = NULL) +
  scale_linetype_manual(values = c("Unadjusted" = "solid", "Adjusted" = "22"),
                        name = NULL) +
  scale_x_continuous("Days") +
  ylab("Cumulative Distribution Function") +
  theme_min(base_family = "sans") +
  theme(legend.position = "none") +
  # labels layer: do NOT inherit global aes (so it doesn't look for adj_f)
  geom_text(
    data = data.table(
      x = c(rep(230,3), rep(30,3)),
      y = c(rep(0.1,3), rep(0.15,3)),
      label = rep(c("IV","PO"), each = 3)
    ),
    aes(x = x, y = y, label = label, colour = label),
    hjust = 1, family = "sans",
    inherit.aes = FALSE
  )+
  guides(linetype = "none", colour = "none", fill = "none") +
  theme(legend.position = "none") 

# -------- RIGHT: PO - IV (with ribbons) --------
rightDT <- plotdt[group == "PO-IV" &
                    analysis %in% c("ITT","per-protocol","as-treated")]

# color to match your PO-IV line color
poiv_col <- "grey"

pRight <- ggplot(rightDT, aes(x = time, y = value, linetype = weighted_f)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = weighted_f),
              alpha = 0.18, colour = NA, show.legend = FALSE) +
  geom_line(colour = "black", linewidth = 1, key_glyph = "path") +                     # center line
  
  # --- add band edges (top and bottom) ---
  geom_line(aes(y = upper), color = poiv_col, linewidth = 0.6) +   # top edge
  geom_line(aes(y = lower), color = poiv_col, linewidth = 0.6) +   # bottom edge
  # ---------------------------------------------------------------

facet_wrap(vars(faclab), ncol = 1) +
  scale_linetype_manual(values = c(Unadjusted = "solid", Adjusted = "22"), name = NULL) +
  scale_fill_manual(values = setNames(rep(poiv_col, 2), c("Unadjusted","Adjusted")),
                    guide = "none") +
  scale_x_continuous(expression(bold("Days"))) +
  scale_y_continuous(position = "right") +
  theme_min(base_family = "sans") +
  theme(axis.title.y = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center")+
  annotate("text", label = "PO \u2212 IV", x = Inf, y = Inf,
           hjust = 1.05, vjust = 1.3)


# -------- Combine and collect legend from right --------
final <- (pLeft | pRight) +
  plot_layout(widths = c(7, 5.4), guides = "collect") &
  theme(legend.position = "bottom",
        legend.justification = "center")

final

if (writeOutput) {
  ggplot2::ggsave(
    filename = "CDFpdfs.pdf",
    plot     = final,
    width    = 12.4,   # inches
    height   = 7,      # inches
    units    = "in",
    device   = grDevices::cairo_pdf
  )
  
}