pt_to_mm = 0.3527777778
mm_to_pt = 1 / pt_to_mm

pubtheme_classic = theme_classic() +
  theme(axis.title = element_text(size = 6,
                                  color = "#000000"), # text size is in points
        axis.text = element_text(size = 5,
                                 color = "#000000"), # text size is in points
        axis.ticks = element_line(color = "#000000",
                                  size = (0.5 * pt_to_mm), # line size is in mm
                                  lineend = "round"),
        axis.line = element_line(color = "#000000",
                                 linewidth = (0.5 * pt_to_mm), # line size in mm
                                 lineend = "round"),
        legend.text = element_text(size = 5, # text size is in points
                                   color = "#000000"),
        legend.title = element_text(size = 6, # text size is in points
                                    color = "#000000"),
        panel.grid = element_line(color = "#FFFFFF"),
        strip.text = element_text(size = 8, color = "#000000")) # text size is in points

pubtheme_bw = theme_bw() +
  theme(axis.title = element_text(size = 6, # text size in points
                                  color = "#000000"), # text size in points
        axis.text = element_text(size = 5, # text size in points
                                 color = "#000000"), # text size in points
        axis.ticks = element_line(color = "#000000",
                                  linewidth = (0.5 * pt_to_mm), # line size in mm
                                  lineend = "round"),
        axis.line = element_line(color = "#000000",
                                 linewidth = (0.5 * pt_to_mm), # line size in mm
                                 lineend = "round"),
        legend.text = element_text(size = 5, # text size in points
                                   color = "#000000"),
        legend.title = element_text(size = 6, # text size in points
                                    color = "#000000"),
        panel.grid = element_line(color = "#FFFFFF"),
        strip.text = element_text(size = 8, color = "#000000"))

pubtheme_cowplot = cowplot::theme_cowplot() +
  theme(axis.title = element_text(size = 6, # text size in points
                                  color = "#000000"), # text size in points
        axis.text = element_text(size = 5, # text size in points
                                 color = "#000000"), # text size in points
        axis.ticks = element_line(color = "#000000",
                                  linewidth = (0.25 * pt_to_mm), # line size in mm
                                  lineend = "round"),
        axis.line = element_line(color = "#000000",
                                 linewidth = (0.5 * pt_to_mm), # line size in mm
                                 lineend = "round"),
        legend.text = element_text(size = 5, # text size in points
                                   color = "#000000"),
        legend.title = element_text(size = 6, # text size in points
                                    color = "#000000"),
        panel.grid = element_line(color = "#FFFFFF"),
        strip.text = element_text(size = 8, color = "#000000"))
