# Load the Volesti library
library(volesti)

A <- rbind(
  c(1, 0, 0),  # x <= 1
  c(-1, 0, 0), # x >= -1
  c(0, 1, 0),  # y <= 1
  c(0, -1, 0), # y >= -1
  c(0, 0, 1),  # z <= 1
  c(0, 0, -1)  # z >= -1
)
b <- rep(1, 6)
polytope <- Hpolytope(A = A, b = b)

samples <- sample_points(polytope, n = 1000, random_walk = "BW")

library(ggplot2)
g <- ggplot(data.frame(x = samples[1, ], y = samples[2, ])) +
  geom_point(aes(x = x, y = y), color = "blue") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  ggtitle("Example Output of Easy Task with Billiard Walk") +
  theme(plot.title = element_text(hjust = 0.5)) # Center the title

print(g)