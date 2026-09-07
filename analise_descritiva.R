grafico_linha <- ggplot(data) + 
  geom_line(aes(x = Data, y = .data[[coluna]]), color = "blue") +
  xlab("Data") + 
  ylab(coluna) +
  ggtitle(paste("Evolução dos Retornos -", coluna))

ggsave(filename = paste0("Grafico_Linha_", coluna, ".png"), plot = grafico_linha, width = 8, height = 4)

# 2. Textos Descritivos no Console
print(paste("================ ANÁLISE DE:", coluna, "================"))
print(summary(retornos))
print(Box.test(retornos, type = "Ljung-Box", lag = 10))

# 3. Gráficos ACF e ACF Quadrado (Lado a Lado)
png(filename = paste0("ACF_", coluna, ".png"), width = 800, height = 400) 

par(mfrow = c(1, 2))
acf(retornos, main = paste("ACF Retornos -", coluna), lag.max = 20)
acf(retornos^2, main = paste("ACF Retornos^2 -", coluna), lag.max = 20)
par(mfrow = c(1, 1))

dev.off() 