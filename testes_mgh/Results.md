# Table comparing results without and with ADOL-C on the MGH test set


| Número do Problema | Nome do Problema             | Número de iterações | Número de avaliações funcionais | Valor da função objetivo | Sup-norma do gradiente projetado | Pontos finais   | Utilizando Adol-c |
|-------------------|------------------------------|---------------------|--------------------------------|--------------------------|---------------------------------|-----------------|------------------|
| 1                 | Rosenbrock           | 61                 | 286                            | 1.300765e-22                    | 7.227963e-11                          | (1.000000e+00, 1.000000e+00)      | Não              |
| 1                | Rosenbrock          | 36                  | 43                            | 5.294921e-21                  | 2.924376e-10                         | (1.000000e+00, 1.000000e+00)      | Sim              |
| 2                 | Freudenstein and Roth           | 46                  | 63                            | 4.898425e+01                   | 8.803225e-07                          | (1.141278e+01, -8.968052e-01)      | Não              |
| 2                 | Freudenstein and Roth           | 46                  | 63                            | 4.898425e+01                   |8.803226e-07                          | (1.141278e+01, -8.968052e-01)      | Sim              |
| 3                 | Powell badly scaled           | 568                  | 2775                            | 2.372760e-07                   |5.833374e-07                          | (1.348146e-05, 7.417596e+00)      | Não              |
| 3                 | Powell badly scaled           | 440                  | 2186                            | 3.398725e-07                  |8.101755e-07                          | (1.375687e-05, 7.269095e+00)      | Sim              |
| 4                 | Brown badly-scaled         | 20                  | 102                            | 8.673617e-19                  |1.862645e-09                          | (1.000000e+06, 2.000000e-06)      | Não              |
| 4                 | Brown badly-scaled         |     20             | 102                            | 2.168404e-19                  |9.313226e-10                          | (1.000000e+06, 2.000000e-06)      | Sim              |
| 14                 | Wood        |       228           | 331                            | 6.142911e-14                  |3.753757e-07                          | (1.000000e+00, 1.000000e+00, 9.999999e-01, 9.999997e-01)      | Não              |
| 14                 | Wood        |       187           | 290                            | 2.529588e-17                  |3.860930e-08                          | (1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00)      | Sim              |
| ...               | ...                          | ...                 | ...                            | ...                      | ...                             | ...             | ...              |
| ...               | ...                          | ...                 | ...                            | ...                      | ...                             | ...             | ...              |
| ...               | ...                          | ...                 | ...                            | ...                      | ...                             | ...             | ...              |
| ...               | ...                          | ...                 | ...                            | ...                      | ...                             | ...             | ...              |
| ...               | ...                          | ...                 | ...                            | ...                      | ...                             | ...             | ...              |
| ...               | ...                          | ...                 | ...                            | ...                      | ...                             | ...             | ...              |
