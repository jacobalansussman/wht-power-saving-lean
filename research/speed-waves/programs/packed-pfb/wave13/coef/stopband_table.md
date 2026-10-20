| proto | case | taps per packed group | quantiser | peak coef | ripple +-dB | worst dB [0.685,1.5) ch | worst dB >= 1.5 ch | median dB >= 1.5 ch | share of stop-band <= -80 dB | coef error floor dB |
|---|---|---|---|---|---|---|---|---|---|---|
| A | packed_s19_L1 | 14 | scale search + active-set minimax search (90 steps, 33 accepted) | 1082 | 0.167 | -60.20 | -75.48 | -92.6 | 98.7 % | -75.5 |
| A | packed_s19_L1 | 7 | scale search + active-set minimax search (90 steps, 52 accepted) | 1477 | 0.166 | -60.27 | -79.36 | -94.1 | 99.7 % | -77.7 |
| A | packed_s19_L1 | 14 | scale search + active-set minimax search (500 steps, 66 accepted) | 1082 | 0.167 | -60.22 | -75.99 | -92.1 | 98.7 % | -75.9 |
| B | packed_s19_L1 | 14 | scale search + active-set minimax search (90 steps, 27 accepted) | 1134 | 0.191 | -61.37 | -75.03 | -92.3 | 98.7 % | -74.6 |
| A | float | - | float (no quantisation) | - | 0.166 | -60.54 | -85.00 | -115.0 | 99.9 % | -318.6 |
| A | w9 | - | round, largest legal scale | 255 | 0.180 | -54.56 | -55.66 | -83.5 | 66.7 % | -55.6 |
| A | w9 | - | scale search (64) | 255 | 0.180 | -54.56 | -55.66 | -83.5 | 66.7 % | -55.6 |
| A | w9 | - | scale search + greedy local search (62 steps) | 255 | 0.179 | -56.37 | -55.91 | -82.4 | 62.9 % | -55.8 |
| A | packed_s18_L1 | 14 | round, largest legal scale | 600 | 0.167 | -60.08 | -67.79 | -89.0 | 90.1 % | -67.7 |
| A | packed_s18_L1 | 14 | scale search (64) | 600 | 0.167 | -60.08 | -67.79 | -89.0 | 90.1 % | -67.7 |
| A | packed_s18_L1 | 14 | scale search + greedy local search (81 steps) | 600 | 0.167 | -60.21 | -69.18 | -87.6 | 89.1 % | -68.7 |
| A | packed_s19_L1 | 14 | round, largest legal scale | 1200 | 0.167 | -59.68 | -74.16 | -93.9 | 99.1 % | -76.5 |
| A | packed_s19_L1 | 14 | scale search (64) | 1081 | 0.167 | -60.26 | -74.19 | -93.0 | 98.7 % | -74.6 |
| A | packed_s19_L1 | 14 | scale search + greedy local search (10 steps) | 1081 | 0.167 | -60.30 | -74.70 | -92.9 | 98.7 % | -75.0 |
| A | packed_s19_L1 | 4 | round, largest legal scale | 1531 | 0.166 | -59.89 | -75.33 | -95.6 | 99.6 % | -76.8 |
| A | packed_s19_L1 | 4 | scale search (64) | 1497 | 0.166 | -60.29 | -76.78 | -95.4 | 99.6 % | -77.7 |
| A | packed_s19_L1 | 4 | scale search + greedy local search (40 steps) | 1497 | 0.166 | -60.19 | -78.44 | -93.7 | 99.7 % | -77.7 |
| A | packed_s19_L1 | 7 | round, largest legal scale | 1482 | 0.166 | -60.17 | -74.91 | -95.2 | 99.5 % | -76.3 |
| A | packed_s19_L1 | 7 | scale search (64) | 1477 | 0.166 | -60.31 | -75.58 | -95.2 | 99.6 % | -77.1 |
| A | packed_s19_L1 | 7 | scale search + greedy local search (27 steps) | 1477 | 0.166 | -60.24 | -77.38 | -93.3 | 99.6 % | -76.7 |
| A | w11 | - | round, largest legal scale | 1023 | 0.167 | -59.78 | -73.17 | -92.7 | 98.1 % | -73.9 |
| A | w11 | - | scale search (64) | 985 | 0.167 | -60.22 | -72.98 | -92.3 | 97.8 % | -72.8 |
| A | w11 | - | scale search + greedy local search (28 steps) | 985 | 0.167 | -60.16 | -73.39 | -91.9 | 97.8 % | -73.2 |
| A | w10 | - | round, largest legal scale | 511 | 0.167 | -59.17 | -65.10 | -87.8 | 86.2 % | -65.1 |
| A | w10 | - | scale search (64) | 508 | 0.167 | -59.65 | -65.30 | -87.9 | 86.1 % | -65.4 |
| A | w10 | - | scale search + greedy local search (80 steps) | 509 | 0.167 | -60.22 | -66.84 | -86.3 | 83.8 % | -66.6 |
| A | w18_round | - | round, peak = full scale | - | 0.166 | -60.53 | -84.94 | -114.9 | 99.9 % | -118.8 |
| A | w12 | - | round, largest legal scale | 2047 | 0.166 | -60.33 | -77.46 | -98.0 | 99.8 % | -79.7 |
| A | w12 | - | scale search (64) | 1956 | 0.166 | -60.19 | -78.99 | -97.5 | 99.8 % | -80.6 |
| A | w12 | - | scale search + greedy local search (19 steps) | 1956 | 0.166 | -60.11 | -79.82 | -97.2 | 99.8 % | -80.4 |
| B | float | - | float (no quantisation) | - | 0.190 | -62.28 | -84.00 | -114.0 | 99.9 % | -318.6 |
| B | packed_s19_L1 | 14 | round, largest legal scale | 1199 | 0.190 | -61.59 | -72.46 | -93.7 | 98.9 % | -74.2 |
| B | packed_s19_L1 | 14 | scale search (64) | 1134 | 0.191 | -61.56 | -74.56 | -93.4 | 98.7 % | -74.3 |
| B | packed_s19_L1 | 14 | scale search + greedy local search (20 steps) | 1134 | 0.191 | -61.45 | -75.19 | -92.3 | 98.7 % | -74.5 |
| B | w18_round | - | round, peak = full scale | - | 0.190 | -62.28 | -84.02 | -114.0 | 99.9 % | -118.1 |
| B | w12 | - | round, largest legal scale | 2047 | 0.190 | -61.70 | -78.03 | -97.8 | 99.8 % | -79.0 |
| B | w12 | - | scale search (64) | 1989 | 0.190 | -61.89 | -78.84 | -97.5 | 99.8 % | -81.1 |
| B | w12 | - | scale search + greedy local search (13 steps) | 1989 | 0.190 | -61.92 | -79.42 | -97.2 | 99.8 % | -80.9 |
