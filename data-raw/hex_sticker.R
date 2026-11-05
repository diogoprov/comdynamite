# Hex sticker for comdynamite with hexSticker. Run from the package root after
# choosing an option (1-4 from hex_options.R, "A"-"E" from hex_options2.R);
# writes man/figures/logo.png. Needs ggplot2, ape and
# hexSticker. The subplots come from data-raw/hex_options.R, drawn in the
# coordinates of the hexagon, so s_x/s_y/s_width/s_height place them where the
# previews show them; adjust p_size if the name does not fit.
library(hexSticker)
option <- 1
FRAME <- FALSE
source("data-raw/hex_options.R")
source("data-raw/hex_options2.R")
looks <- list(
  list(p = g1, fill = "#1D3557", border = "#E63946", text = "#F1FAEE", p_y = 0.5),
  list(p = g2, fill = "#264653", border = "#E9C46A", text = "#E9C46A", p_y = 0.47),
  list(p = g3, fill = "#FFB703", border = "#023047", text = "#023047", p_y = 0.47),
  list(p = g4, fill = "#F1FAEE", border = "#1D3557", text = "#1D3557", p_y = 0.5),
  A = list(p = pA, fill = "#264653", border = "#E9C46A", text = "#E9C46A", p_y = 0.5),
  B = list(p = pB, fill = "#F1FAEE", border = "#1D3557", text = "#1D3557", p_y = 0.52),
  C = list(p = pC, fill = "#1D3557", border = "#F4A261", text = "#F4A261", p_y = 0.5),
  D = list(p = pD, fill = "#2A9D8F", border = "#264653", text = "#F1FAEE", p_y = 0.5),
  E = list(p = pE, fill = "#264653", border = "#E9C46A", text = "#E9C46A", p_y = 0.47)
)
look <- looks[[option]]
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
sticker(look$p + theme_transparent(), package = "comdynamite",
        p_size = 17, p_color = look$text, p_y = look$p_y, p_fontface = "bold",
        s_x = 1, s_y = 1, s_width = 1.76, s_height = 2.04,
        h_fill = look$fill, h_color = look$border, h_size = 1.6,
        filename = "man/figures/logo.png", dpi = 300)
