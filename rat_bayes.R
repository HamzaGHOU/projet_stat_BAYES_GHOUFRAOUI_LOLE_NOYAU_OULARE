

###on charge les data 
N <-30
T <- 5
y <-structure(c(151, 145, 147, 155, 135, 159, 141, 159, 177, 134, 
              160, 143, 154, 171, 163, 160, 142, 156, 157, 152, 154, 139, 146, 
              157, 132, 160, 169, 157, 137, 153, 199, 199, 214, 200, 188, 210, 
              189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 212, 
              203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, 
              263, 237, 230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, 
              248, 234, 243, 259, 246, 253, 225, 229, 250, 237, 257, 261, 248, 
              219, 244, 283, 293, 312, 272, 280, 298, 275, 297, 350, 260, 313, 
              273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, 
              286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 
              338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 
              334, 302, 302, 323, 331, 345, 333, 316, 291, 324), .Dim = c(30, 
                                                                          5))
x <- c(8.0, 15.0, 22.0, 29.0, 36.0)
xbar <- 22.0
####init 
alpha <- c(250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 
           250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250)
beta <- c(6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)
mu_alpha <- 150
mu_beta <- 10
sigmasq_y <- 1
sigmasq_alpha <- 1
sigmasq_beta <- 1

##graph intro 

library(ggplot2)
library(tidyr)
data_long <- data.frame(x = rep(x, nrow(y)), y = c(t(y)), Rat = factor(rep(1:nrow(y), each = length(x))))

palette_couleurs <- rainbow(length(unique(data_long$Rat)))

ggplot(data_long, aes(x = x, y = y, group = Rat, color = Rat, linetype = Rat)) +
  geom_path() +
  geom_point(shape = 1) +
  labs(x = "Jour X_j", y = "Poids") 
scale_color_manual(values = palette_couleurs) +  # Utiliser une palette de couleurs
  scale_linetype_manual(values = rep("dashed", length(unique(data_long$Rat))))  # Utiliser des lignes pointillées

### ech de gibbs 

rat_gibbs <- function(xbar,       ,y,x)