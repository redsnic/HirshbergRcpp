
library(dplyr)
library(tidyr)
library(ggplot2)

Ns <- c(1000,3000,5000,7000, 10000, 20000, 30000, 50000, 70000, 100000)

tara <- 75964
NWSW_mem <- c(96956, 100976, 193648, 245284, 1593380, NA, NA, NA, NA, NA)
HIRSH_mem <- c(75964, 75964, 75964, 75964, 75964, 75964, 75964, 75964, 82072, 89260)

df_mem <- data.frame(Ns, NWSW_mem-tara, HIRSH_mem-tara)
colnames(df_mem) <- c("N", "NWSW", "Hirshberg")
df_mem <- df_mem %>% gather("NWSW","Hirshberg", key="algorithm", value="KBytes")

space.plot <- (
ggplot(df_mem) +
  geom_line( aes(x = N, y=KBytes, color=algorithm)) +
  scale_x_log10() +
  labs(title = "Complessità spaziale", x="Dimensione dell'istanza", y="Utilizzo in KB della memoria")
)

load("benchmark.Rds")

total_sizes <- c()
total_modes <- c()
total_algorithms <- c()
total_values <- c()

for(i in 1:length(Ns)){
  
  thg <- benchs_global[[i]]$th
  tnwg <- benchs_global[[i]]$tnw
  thl <- benchs_local[[i]]$th
  tnwl <- benchs_local[[i]]$tnw
  
  problem_sizes <- rep(Ns[i], times = length(thg)+length(tnwg)+length(thl)+length(tnwl))
   
  values <- c(thg, tnwg, thl, tnwl)
  modes <- c(
    rep("global", times = length(thg)+length(tnwg)),
    rep("local", times = length(thl)+length(tnwl))
  )
  algorithms <- c(
    rep("Hirshberg", times = length(thg)),
    rep("NWSW", times = length(tnwg)),
    rep("Hirshberg", times = length(thl)),
    rep("NWSW", times = length(tnwl))
  )
  
  total_sizes <- c(total_sizes, problem_sizes)
  total_modes <- c(total_modes, modes)
  total_algorithms <- c(total_algorithms, algorithms)
  total_values <- c(total_values,values)
  
}
 
df_time <- data.frame(total_sizes, total_modes, total_algorithms, total_values)
colnames(df_time) <- c("dimensione", "mode", "algoritmo", "tempo")

time.plot <- (
ggplot(df_time) +
  geom_boxplot( aes(x = as.factor(dimensione), y=tempo, linetype = algoritmo, fill = mode ) ) + 
  labs(title = "Complessità temporale", x="Dimensione dell'istanza", y="tempo (ms)")
)

time.plot.log <- (
  ggplot(df_time) +
    geom_boxplot( aes(x = as.factor(dimensione), y=tempo, linetype = algoritmo, fill = mode ) ) + 
    labs(title = "Complessità temporale", x="Dimensione dell'istanza", y="tempo (ms)") + 
    scale_y_log10()
)



