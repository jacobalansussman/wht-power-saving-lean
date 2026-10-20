| file | test | degenerate | trivial | excess | non-trivial collapse |
|---|---|---|---|---|---|
| out_annih1_a.txt | annih_camp menu=xw(55) list_annih1.txt[0 | ? | ? | ? | - |
| out_annih1_a.txt | annih_camp chain=1;1;2 m=4 | 2 | 2 | 0 | - |
| out_annih1_a.txt | annih_camp chain=1;1;2 m=7 | 0 | 0 | 0 | - |
| out_annih1_a.txt | annih_camp chain=1;1;2 m=5 | 0 | 0 | 0 | - |
| out_annih1_a.txt | annih_camp chain=1;1;2 m=6 | 1 | 1 | 0 | - |
| out_annih1_b.txt | annih_camp menu=xw(55) list_annih1.txt[50 | ? | ? | ? | - |
| out_annih1_b.txt | annih_camp chain=1;a0;2 m=4 | 0 | 0 | 0 | - |
| out_annih1_b.txt | annih_camp chain=1;1;2,a0 m=4 | 2 | 2 | 0 | - |
| out_annih1_b.txt | annih_camp chain=1;1;2 m=4,a0 | 2 | 2 | 0 | - |
| out_annih1_b.txt | annih_camp chain=1,a0i;1;2 m=4 | 2 | 2 | 0 | - |
| out_annih1_c.txt | annih_camp menu=xw(55) list_annih1.txt[100 | ? | ? | ? | - |
| out_annih1_c.txt | annih_camp chain=1;1;2 m=4,A.T2 | 0 | 0 | 0 | - |
| out_annih1_c.txt | annih_camp chain=1,T2.A;1;2 m=4 | 2 | 2 | 0 | - |
| out_annih1_c.txt | annih_camp chain=1;T2.A;2 m=4 | 0 | 0 | 0 | - |
| out_annih1_c.txt | annih_camp chain=1;1;2,T2.A m=4 | 0 | 0 | 0 | - |
| out_annih1_c.txt | annih_camp chain=1;1;2 m=4,T2.A | 0 | 0 | 0 | - |
| out_annih1_d.txt | annih_camp menu=xw(55) list_annih1.txt[150 | ? | ? | ? | - |
| out_annih1_d.txt | annih_camp chain=1;T2.a2i;2 m=4 | 0 | 0 | 0 | - |
| out_annih1_d.txt | annih_camp chain=1;1;2,T2.a2i m=4 | 0 | 0 | 0 | - |
| out_annih1_d.txt | annih_camp chain=1;1;2 m=4,T2.a2i | 0 | 0 | 0 | - |
| out_annih1_d.txt | annih_camp chain=1,a0.T2.a0i;1;2 m=4 | 2 | 2 | 0 | - |
| out_annih3_longC.txt | annih3 moves=4;4,T2 R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 2070 rows 4608 | 2 | 2 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,T2 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,A.T2.Ai R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 2070 rows 4608 | 2 | 2 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,A.T2.Ai R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,CS01 R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 2070 rows 4608 | 2 | 2 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,CS01 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 co | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,Ai.CS01.A R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 2070 rows 46 | 2 | 2 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,Ai.CS01.A R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,lT2 R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 2070 rows 4608 | 2 | 2 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,lT2 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 col | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;7 R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 2070 rows 4608 | 2 | 2 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;7 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 23 | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4,Pr0;4,A R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 cols 2070 rows 4608 | 2 | 2 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4,Pr0;4,A R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2 FS=[3, 4] S=1 c | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,T2 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1,T2;1;2 FS=[3, 4] S=1 c | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,T2 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;A.T2.Ai;2 FS=[3, 4] S= | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,T2 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1;1;2,A.CS01.Ai FS=[3, 4 | 0 | 0 | 0 | 0 |
| out_annih3_longC.txt | annih3 moves=4;4,T2 R=['1', 'T2', 'CS01', 'A.T2.Ai', 'Ai.CS01.A', 'Pr0', 'A', 'Ai', 'lT2', 'E01'] C=xw+(115) B=xw+(115) chain=1,A;Pr0;2 FS=[3, 4] S=1  | 0 | 0 | 0 | 0 |
| out_annih5_batch.txt | annih4 moves=1;2;4|1 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 14 | 14 | 0 | 1 |
| out_annih5_batch.txt | annih4 moves=1;2;4|4,T2 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 1 | 1 | 0 | 0 |
| out_annih5_batch.txt | annih4 moves=1;2;4|7 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 1 | 1 | 0 | 0 |
| out_annih5_batch.txt | annih4 moves=1;2,T2;4|4,A.T2.Ai |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 0 | 0 | 0 | 0 |
| out_annih5_batch.txt | annih4 moves=1;2;4,CS01|4,Ai.CS01.A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 1 | 1 | 0 | 0 |
| out_annih5_batch.txt | annih4 moves=4;2;1|1,A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 1 | 1 | 0 | 0 |
| out_annih5b.txt | annih4 moves=1;2;4|4,T0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 1 | 1 | 0 | 0 |
| out_annih5b.txt | annih4 moves=1;2;4|4,CS12 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 1 | 1 | 0 | 0 |
| out_annih5b.txt | annih4 moves=1;2;4|7,A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 1 | 1 | 0 | 0 |
| out_annih5b.txt | annih4 moves=1,T2;2,T0;4,CS12|4,A.T0.Ai |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2 | 0 | 0 | 0 | 0 |
| out_annih7.txt | annih3 moves=4;4,T0;4,T1;4,CS12;4,CS02;4,CCZ R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw(55) B=xw(55) chain=1;1;2 FS=[3, 4] S=2 col | 2 | 2 | 0 | 0 |
| out_annih7.txt | annih3 moves=4;4,A.T0.Ai;4,Ai.CS12.A;4,lT0;4,Pr2;4,B.T0.Bi R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw(55) B=xw(55) chain=1;1;2 FS= | 2 | 2 | 0 | 0 |
| out_annih7.txt | annih3 moves=4;4,T0;4,CS02;7;4,A.CS02.Ai;4,Ai.T1.A R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw(55) B=xw(55) chain=1,T0;1;2,CS12 FS= | 0 | 0 | 0 | 0 |
| out_annih1_m2_4.txt | annih_camp menu=xw(55) list_annih1.txt[4 | ? | ? | ? | - |
| out_annih1_m2_4.txt | annih_camp chain=1;1;2 m=3 | 0 | 0 | 0 | - |
| out_annih1_m2_4.txt | annih_camp chain=1;1;3 m=4 | 0 | 0 | 0 | - |
| out_annih1_mac_25.txt | annih_camp menu=xw(55) list_annih1.txt[25 | ? | ? | ? | - |
| out_annih1_mac_25.txt | annih_camp chain=1,CS12;1;2 m=4 | 2 | 2 | 0 | - |
| out_annih1_mac_25.txt | annih_camp chain=1;CS12;2 m=4 | 0 | 0 | 0 | - |
| out_annih1_mac_25.txt | annih_camp chain=1;1;2,CS12 m=4 | 0 | 0 | 0 | - |
| out_annih1_mac_25.txt | annih_camp chain=1;1;2 m=4,CS12 | 0 | 0 | 0 | - |
| out_annih1_mac_25.txt | annih_camp chain=1,CS02;1;2 m=4 | 2 | 2 | 0 | - |
| out_annih1_mac_25.txt | annih_camp chain=1;CS02;2 m=4 | 0 | 0 | 0 | - |
| out_annih1_mac_25.txt | annih_camp chain=1;1;2,CS02 m=4 | 2 | 2 | 0 | - |
| out_annih1_mac_25.txt | annih_camp chain=1;1;2 m=4,CS02 | 0 | 0 | 0 | - |
| out_annih3_longord.txt | annih3 moves=4 R=['1', 'A.T2.B', 'A.B', 'A.C01.B', 'a0.T2.b0', 'A.T0.B', 'A.CS01.B', 'B.T2.A', 'a0.b0', 'A.T2.B.T0'] C=xw(55) B=xwg(82) chain=1;1;2 FS | 0 | 0 | 0 | - |
| out_annih3_multi.txt | annih3 moves=4;7;5;6;3;4,T2;4,A;4,A.T2.Ai;4,CS01 R=['1', 'X1', 'X2', 'X4', 'Z1', 'Z2', 'Z4', 'P01', 'P12'] C=xw(55) B=xw(55) chain=1;1;2 FS=[3, 4] S=2 | 2 | 2 | 0 | - |
| out_annih4.txt | annih4 moves=1;2;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 3 | 3 | 0 | - |
| out_annih4_batch.txt | annih4 moves=1;2;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 3 | 3 | 0 | - |
| out_annih4_batch.txt | annih4 moves=1;2;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih4_batch.txt | annih4 moves=1;2,T2;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1;2,T2;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1;2;4,A.T2.Ai |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 3 | 3 | 0 | 1 |
| out_annih4_batch.txt | annih4 moves=1;2;4,A.T2.Ai |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1,A;2,CS01;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1,A;2,CS01;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1;2;7 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 2 | 2 | 0 | 1 |
| out_annih4_batch.txt | annih4 moves=1;2;7 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1,T2;2,A;4,CS01 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1,T2;2,A;4,CS01 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1;2,Ai.CS01.A;4,T0 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1;2,Ai.CS01.A;4,T0 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih4_batch.txt | annih4 moves=1,B;2;4,A |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 1 | 0 | 1 | 1 |
| out_annih4_batch.txt | annih4 moves=1,B;2;4,A |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | 0 |
| out_annih6.txt | annih6 chain=1;2 t0=4 t1=4 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 24 | 24 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=4 t1=7 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=4 t1=7 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=4 t1=4,T2 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 24 | 24 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=4 t1=4,T2 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 23 | 23 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=4,A.T2.Ai t1=4 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 24 | 24 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=4,A.T2.Ai t1=4 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 23 | 23 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=4 t1=1 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 2 | 2 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=4 t1=1 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 1 | 1 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=4 t1=2,A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=4 t1=2,A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=7 t1=3 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=7 t1=3 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=4,CS01 t1=4,Ai.CS01.A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 24 | 24 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=4,CS01 t1=4,Ai.CS01.A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 23 | 23 | 0 | - |
| out_annih6.txt | annih6 chain=1;2 t0=4 t1=6,T0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6.txt | annih6 chain=1,T2;2,A t0=4 t1=6,T0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1;2 t0=4 t1=4,T0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1,T0;2,CS02 t0=4 t1=4,T0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1;2 t0=4,T0 t1=4,CS12 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1,T0;2,CS02 t0=4,T0 t1=4,CS12 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1;2 t0=4,A.T0.Ai t1=4,Ai.CS12.A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1,T0;2,CS02 t0=4,A.T0.Ai t1=4,Ai.CS12.A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1;2 t0=4,CCZ t1=7 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1,T0;2,CS02 t0=4,CCZ t1=7 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1;2 t0=4,T1 t1=4,B.T0.Bi |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1,T0;2,CS02 t0=4,T1 t1=4,B.T0.Bi |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1;2 t0=4,Pr2 t1=4,lT0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6b.txt | annih6 chain=1,T0;2,CS02 t0=4,Pr2 t1=4,lT0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,CS01;2,CS12 t0=4 t1=4,T0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,T2;2,T0 t0=4 t1=4,T0 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,CS01;2,CS12 t0=4,T0 t1=4,CS12 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,T2;2,T0 t0=4,T0 t1=4,CS12 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,CS01;2,CS12 t0=4,A.T0.Ai t1=4,Ai.CS12.A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,T2;2,T0 t0=4,A.T0.Ai t1=4,Ai.CS12.A |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,CS01;2,CS12 t0=4,T1 t1=4,B.T0.Bi |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih6c.txt | annih6 chain=1,T2;2,T0 t0=4,T1 t1=4,B.T0.Bi |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 0 | 0 | 0 | - |
| out_annih_bigB.txt | annih2 R=xw(55) C=xw(55) B=xw+(115) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 2119, rank A3 mod L2 2117 | 2 | 2 | 0 | - |
| out_annih_first.txt | annih2 R=wd(28) C=wd(28) B=wd(28) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 609, rank A3 mod L2 607 | 2 | ? | ? | - |
| out_annih_first.txt | annih2 R=xqm(55) C=xqm(55) B=xqm(55) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 1818, rank A3 mod L2 1816 | 2 | ? | ? | - |
| out_annih_first.txt | annih2 R=xw(55) C=xw(55) B=xw(55) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 2119, rank A3 mod L2 2117 | 2 | ? | ? | - |
| out_annih_lvl0.txt | annih2 R=eng(22) C=eng(22) B=eng(22) chain=id m=(1,1) FS=[3, 4] rows=4608 rank A3 387, rank A3 mod L2 387 | 0 | 0 | 0 | - |
| out_annih_lvl0.txt | annih2 R=xw(55) C=xw(55) B=xw(55) chain=id m=(1,1) FS=[3, 4] rows=4608 rank A3 2445, rank A3 mod L2 2445 | 0 | 0 | 0 | - |
| out_annih_lvl1.txt | annih2 R=eng(22) C=eng(22) B=eng(22) chain=1 m=(2,1) FS=[3, 4] rows=4608 rank A3 405, rank A3 mod L2 403 | 2 | 2 | 0 | - |
| out_annih_lvl1.txt | annih2 R=xw(55) C=xw(55) B=xw(55) chain=1 m=(2,1) FS=[3, 4] rows=4608 rank A3 2513, rank A3 mod L2 2511 | 2 | 2 | 0 | - |
| out_annih_lvl1.txt | annih2 R=eng(22) C=eng(22) B=eng(22) chain=1 m=(1,1) FS=[3, 4] rows=4608 rank A3 387, rank A3 mod L2 365 | 22 | 22 | 0 | - |
| out_annih_lvl1.txt | annih2 R=xw(55) C=xw(55) B=xw(55) chain=1 m=(1,1) FS=[3, 4] rows=4608 rank A3 2445, rank A3 mod L2 2393 | 52 | 52 | 0 | - |
| out_annih_orbit.txt | annih3 moves=4 R=powers of A.B (0..7) C=xw(55) B=xw(55) chain=1;1;2 FS=[3, 4] S=1 cols 440 rows 4608 | 0 | 0 | 0 | 0 |
| out_annih_orbit.txt | annih3 moves=4 R=powers of A.T2.B (0..7) C=xw(55) B=xw(55) chain=1;1;2 FS=[3, 4] S=1 cols 440 rows 4608 | 0 | 0 | 0 | 0 |
| out_nt.txt | annih2 R=eng(22) C=eng(22) B=eng(22) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 360, rank A3 mod L2 358 | 2 | 2 | 0 | 0 |
| out_nt.txt | annih2 R=wd(28) C=wd(28) B=wd(28) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 609, rank A3 mod L2 607 | 2 | 2 | 0 | 0 |
| out_nt.txt | annih2 R=xqm(55) C=xqm(55) B=xqm(55) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 1818, rank A3 mod L2 1816 | 2 | 2 | 0 | 0 |
| out_nt.txt | annih2 R=xw(55) C=xw(55) B=xw(55) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 2119, rank A3 mod L2 2117 | 2 | 2 | 0 | 0 |
| out_nt.txt | annih2 R=xw(55) C=xw(55) B=xw+(115) chain=1;1;2 m=(4,1) FS=[3, 4] rows=4608 rank A3 2119, rank A3 mod L2 2117 | 2 | 2 | 0 | 0 |
| out_nt.txt | annih2 R=xw(55) C=xw(55) B=xw(55) chain=1 m=(2,1) FS=[3, 4] rows=4608 rank A3 2513, rank A3 mod L2 2511 | 2 | 2 | 0 | 0 |
| out_nt.txt | annih2 R=xw(55) C=xw(55) B=xw(55) chain=1 m=(1,1) FS=[3, 4] rows=4608 rank A3 2445, rank A3 mod L2 2393 | 52 | 52 | 0 | 4 |
| out_nt.txt | annih4 moves=1;2;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1 | 3 | 3 | 0 | 1 |
| out_nt.txt | annih6 chain=1;2 t0=4 t1=4 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 24 | 24 | 0 | 4 |
| out_nt.txt | annih6 chain=1;2 t0=4 t1=1 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1 | 2 | 2 | 0 | 0 |

142 tests; with EXCESS > 0: 1; with non-trivial collapse > 0: 7
   ('out_annih4_batch.txt', 'annih4 moves=1,B;2;4,A |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1', '1', '0', '1', '1')
   ('out_annih5_batch.txt', 'annih4 moves=1;2;4|1 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=2', '14', '14', '0', '1')
   ('out_annih4_batch.txt', 'annih4 moves=1;2;4,A.T2.Ai |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1', '3', '3', '0', '1')
   ('out_annih4_batch.txt', 'annih4 moves=1;2;7 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1', '2', '2', '0', '1')
   ('out_annih4_batch.txt', 'annih4 moves=1,B;2;4,A |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1', '1', '0', '1', '1')
   ('out_nt.txt', 'annih2 R=xw(55) C=xw(55) B=xw(55) chain=1 m=(1,1) FS=[3, 4] rows=4608 rank A3 2445, rank A3 mod L2 2393', '52', '52', '0', '4')
   ('out_nt.txt', 'annih4 moves=1;2;4 |R|=9 |H|=30 |K|=9 |B|=30 FS=[3, 4] S=1', '3', '3', '0', '1')
   ('out_nt.txt', 'annih6 chain=1;2 t0=4 t1=4 |R|=9 |H|=30 |K|=5 |B|=30 FS=[3, 4] S=1', '24', '24', '0', '4')
