# creating figure 2 and figure S2
library(patchwork)

source("inference_comparison_fw5.R")
source("inference_comparison_switching.R")
source("inference_comparison_logistic.R")

fig2 <- plot_fig2a + plot_fig2b + plot_fig2c + plot_fig2d + plot_fig2e + plot_fig2f + 
  plot_fig2g + plot_fig2h + plot_fig2i + plot_fig2j + plot_fig2k + plot_fig2l + 
  plot_layout(ncol = 3)
ggsave("fig/figure2.png", dpi = 300, width = 11, height = 10)

figS2 <- plot_figS2a + plot_figS2b + plot_figS2c + plot_figS2d + plot_figS2e + plot_figS2f + 
  plot_figS2g + plot_figS2h + plot_figS2i + plot_figS2j + plot_figS2k + plot_figS2l + 
  plot_layout(ncol = 3)
ggsave("fig/figureS2.png", dpi = 300, width = 11, height = 10)
