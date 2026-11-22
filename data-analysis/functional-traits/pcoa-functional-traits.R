##PCoA (numeric and categorical variables)
library(GGally)
library(performance)
library(vegan)
library(patchwork)
library(DHARMa)
library(glmmTMB)
library(car)
#data 
setwd('~/projects/data-analysis/functional-traits/data/')

pcoa_traits <- read.csv('pca_traits.csv', header = T, sep = ',', stringsAsFactors = T, fileEncoding = 'latin1')

#CHEQUEO DE COLINEALIDAD
variable_ficticea <- rep(1, 22)
matrix_col <- cbind(variable_ficticea, pcoa_traits[,c(9:18)])
modelo <- lm(variable_ficticea ~  fenologia_hoja_ID + conc_P + AFefc + AFE + FMF + AF_AC + Ft + CMSh + Dl + CAS,  data = matrix_col)
colinealidad <- check_collinearity(modelo)
colinealidad

#elimino CAS que tiene vif alto
modelo <- lm(variable_ficticea ~  fenologia_hoja_ID + conc_P + AFefc + AFE + FMF + AF_AC + Ft + CMSh + Dl,  data = matrix_col)
colinealidad <- check_collinearity(modelo)
colinealidad

##elimino AF/AC que tiene vif alto
modelo <- lm(variable_ficticea ~  fenologia_hoja_ID + conc_P + AFefc + AFE + FMF + Ft + CMSh + Dl,  data = matrix_col)
colinealidad <- check_collinearity(modelo)
colinealidad

traits <- pcoa_traits[, c(9:17)] 
head(traits)
colnames(traits)

dist<-vegdist(traits, method="gower", binary=F) #c?lculo de distancia de Gower hace estandarizacion de variables cuantitativas

PCoA<-cmdscale(dist, eig=T, add=T)
PCoA

## calculo de la variablidad que absorbe cada uno de los ejes 
PCoA$eig/sum(PCoA$eig)*100

# calculo de la importancia relativa de las variables en la ordenacion
which(dist<0) #para saber si hay valores negativos. Si hubiese hay que hacer correcci?n
which(traits<0)

vectores<-envfit(PCoA, traits, perm=1000)
vectores
rownames(vectores$vectors$arrows) <- c('Fenología foliar', 'P', 'AFefc', 'AFE',
                                       'FMF', 'AF/AC', 'Ft', 'CMSh', 'Dl')

# GR?FICO
windows()
## Mapa de colores para los tipos de dormici?n
color_mapping <- c("Invasora" = "#717171", "No invasora" = "#D0D0D0")


pcoa_coords <- as.data.frame(PCoA$points[, 1:2]) # Primeros dos ejes
colnames(pcoa_coords) <- c("PCo1", "PCo2")
pcoa_coords$status <- pcoa_traits$status # Invasoras/no invasoras
pcoa_coords$par<- pcoa_traits$par

# Como los puntos 5 (no inv) y 6 (inv) se superponen
pcoa_coords$PCo1[6] <- pcoa_coords$PCo1[6] + 0.03 # Desplazar ligeramente en el eje X
pcoa_coords$PCo2[6] <- pcoa_coords$PCo2[6] - 0.03  # Desplazar ligeramente en el eje Y

#tiff(filename = "PCoA_without_legend.tiff",
#     width = 17.4, height = 14, units = "cm", res = 300)

png(filename = "PCoA1_without_legend.png",
    width = 17.4, height = 14, units = "cm", res = 300)

par(mar = c(4, 4, 2, 2))
fig <- ordiplot(PCoA, choices = c(1, 2), type = "none", display = "sites",  
                xlab = "PCo 1 (26.36 %)", ylab = "PCo 2 (15.25 %)")
plot.window(xlim = c(-0.7, 0.7), ylim = c(-0.5, 0.5))
abline(h = 0, col = "grey70", lty = 2)
abline(v = 0, col = "grey70", lty = 2)
# Diagonales
#segments(-0.7, 0.5, 0.7, -0.5, col = "grey80", lty = 2)
#segments(-0.7, -0.5, 0.7, 0.5, col = "grey80", lty = 2)
points(pcoa_coords$PCo1, pcoa_coords$PCo2,
       pch = 16,
       col = color_mapping[pcoa_coords$status], # Color seg?n status
       cex = 1)
plot(vectores, p.max = 0.05, col="black", cex=0.9, font=1)
text(pcoa_coords$PCo1, pcoa_coords$PCo2, 
     labels = pcoa_coords$par, 
     col = "black",
     font = 1,
     cex = 0.9, 
     pos = 2,
     offset = 0.5)
dev.off()



#MODELOS
axis1 <- PCoA$points[,1]
axis2 <- PCoA$points[,2]
pcoa_scores <- data.frame(axis1, axis2, pcoa_traits$status, pcoa_traits$especie, pcoa_traits$par)
pcoa_scores$pcoa_traits.especie <- as.factor(pcoa_scores$pcoa_traits.especie)
pcoa_scores$pcoa_traits.par <- as.factor(pcoa_scores$pcoa_traits.par)
colnames(pcoa_scores) <- c('axis1', 'axis2', 'status', 'especie', 'par')
str(pcoa_scores)

fit_pcoa1 <- glmmTMB(axis1 ~ status + (1|status/especie) + (1|par/especie), data = pcoa_scores)
simulationOutput <- simulateResiduals(fittedModel = fit_pcoa1, plot = F)
plot(simulationOutput)
Anova(fit_pcoa1, type = 'II')

fit_pcoa2 <- glmmTMB(axis2 ~ status + (1|status/especie) + (1|par/especie) , data = pcoa_scores)
simulationOutput <- simulateResiduals(fittedModel = fit_pcoa2, plot = F)
plot(simulationOutput)
Anova(fit_pcoa2, type = 'II')

#boxplots de las spp seg?n ststus en la puntuaci?n de PCo1 y PCo 2
plot_pco1 <- ggplot(pcoa_scores, aes(x = status, y = axis1, fill = status)) + 
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  labs(y = 'PCo 1', x = 'status', tag = 'a') +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 2, shape = 19) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  theme_bw() + 
  theme(legend.position = 'none')
plot_pco1


df <- data.frame(status = 'Invasora',
                 y = 0.30, 
                 text = '***')
plot_pco2 <- ggplot(pcoa_scores, aes(x = status, y = axis2, fill = status)) + 
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  labs(y = 'PCo 2', x = 'status', tag = 'b') +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 2, shape = 19) +
  geom_text(data = df, aes(x = status, y = y, label = text), size = 7, position = position_nudge(x = 0.5)) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  theme_bw() +
  theme(legend.position = 'none', plot.margin = unit(c(0, 0, 0, 1), "cm"))
plot_pco2

plots <- plot_pco1 + plot_pco2
plots
