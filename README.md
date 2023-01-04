# Masters Thesis

Reproducible code for my Masters thesis (Spring 2023)

# Simulation

| Parameter | Explanation                                  | Choices                                  |
|-----------|----------------------------------------------|------------------------------------------|
| f(x, y)   | Bivariate generating distribution            | {Normal, Ordinal}                        |
| ρ         | Correlation                                  | {-0.9,   -0.5, -0.25, 0, 0.25, 0.5, 0.9} |
| δ         | True difference in means                     | {0, 0.25, 0.5}                           |
| n         | Sample size (per group)                      | {10, 20, 50, 100, 200}                   |
| m/n       | Proportion of matched samples                | {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5}      |
| σ_x, σ_y  | Variance                                     | {1}                                      |
| R         | Number of repetitions per simulated scenario | 10000                                    |