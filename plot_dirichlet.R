plot_dirichlet <- function(data_matrix){
  library(ggtern)
  ggtern(data = data.frame(x = data_matrix[,1], y = data_matrix[,2], z = data_matrix[,3]), aes(x, y, z))+
    stat_density_tern(geom = 'polygon', n = 500, 
                      aes(fill  = after_stat(level), alpha = after_stat(level)),
                      bdl = 0.01) + 
    scale_fill_gradient(low = "blue", high = "red", name = "", breaks = 1:5, 
                        labels = c("low", "", "", "", "high"))  +
    scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
    scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
    scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
    guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
    theme_rgbg() +
    theme_noarrows() 
}