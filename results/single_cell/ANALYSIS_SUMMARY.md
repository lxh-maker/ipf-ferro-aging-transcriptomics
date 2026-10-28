# GSE135893 targeted single-cell analysis summary

- Strict cohort: 12 IPF donors, 10 controls.
- Cells retained: 89,326 (IPF and Control only).
- Author-provided cell types: 31.
- Ferro-aging genes used in targeted composite score: 95.
- Differential tests are donor-level pseudobulk tests; cells are not treated as independent biological replicates.
- The targeted ferro-aging score is a standardized composite over the supplied gene set, not a full-transcriptome AUCell/UCell score.

## Five-gene compartment results

```text
  gene  p_value  median_IPF  median_Control      FDR
 PDE4B 0.766679    5.628934        4.915081 0.803316
  LIFR 0.803316    4.408921        4.169794 0.803316
S100A8 0.575157    9.497929        9.083556 0.803316
  EGR1 0.276605    7.854656        8.307130 0.691512
  LCN2 0.009199    8.113289        5.282701 0.045997
```

## Top candidate/mechanism pseudobulk results

```text
          celltype   gene  n_IPF_samples  n_Control_samples  median_IPF_log2CPM  median_Control_log2CPM  IPF_minus_Control  p_value      FDR
          SCGB3A2+    LIF             12                 10            2.161916                0.000000           2.161916 0.000123 0.087346
               AT2  HIF1A             12                 10            6.731997                5.466050           1.265946 0.000326 0.087346
               AT2  ACSL4             12                 10            7.587340                5.716553           1.870787 0.001384 0.100858
         Monocytes PLA2G7             12                 10            5.932607                3.690358           2.242248 0.001715 0.100858
 Endothelial Cells   LIFR             12                  9            8.434844                9.621785          -1.186942 0.001992 0.108070
  Transitional AT2  HIF1A             11                 10            5.890627                4.317120           1.573506 0.002592 0.121502
       Macrophages PLA2G7             12                 10            7.235715                5.005099           2.230616 0.003344 0.126185
    Myofibroblasts  HIF1A             12                  9            7.905622                6.869961           1.035661 0.003999 0.132672
       Macrophages S100A8             12                 10            9.294616                7.810363           1.484253 0.007574 0.172449
 Endothelial Cells S100A9             12                  9            4.153371                6.332631          -2.179260 0.009488 0.187177
       Macrophages S100A9             12                 10            8.968147                7.770842           1.197305 0.013411 0.218598
            MUC5B+  PTGS2             12                  9            3.164464                0.000000           3.164464 0.014517 0.221066
          SCGB3A2+   LIFR             12                 10            4.804863                5.694159          -0.889296 0.014671 0.221470
  Transitional AT2  ACSL4             11                 10            5.576925                3.683062           1.893864 0.016305 0.236795
PLIN2+ Fibroblasts  ACSL4              8                  3            8.275558                0.000000           8.275558 0.018634 0.240323
  Transitional AT2 S100A8             11                 10            0.000000                2.322409          -2.322409 0.019228 0.240323
          SCGB3A2+  SPHK1             12                 10            0.000000                0.000000           0.000000 0.027812 0.284583
          SCGB3A2+  PDE4B             12                 10            1.452883                0.000000           1.452883 0.029085 0.289430
       Fibroblasts  PDE4B             11                  9            5.243415                0.000000           5.243415 0.029742 0.291201
       Macrophages  SPHK1             12                 10            5.024992                2.801086           2.223906 0.032114 0.296544
```

## Top ferro-aging score shifts

```text
                celltype  n_IPF_samples  n_Control_samples  median_IPF  median_Control  IPF_minus_Control  p_value      FDR
                    cDCs             12                 10   -0.025529        0.062762          -0.088291 0.051754 0.580323
Differentiating Ciliated             12                  7   -0.083655       -0.156706           0.073052 0.068350 0.580323
                     AT1              8                  9    0.043995        0.006242           0.037753 0.074455 0.580323
   Proliferating T Cells             10                  8    0.014729       -0.097482           0.112211 0.101102 0.580323
                Ciliated             12                 10   -0.031068       -0.058641           0.027574 0.137913 0.580323
                  MUC5B+             12                  9   -0.120456       -0.159772           0.039316 0.165807 0.580323
       SCGB3A2+ SCGB1A1+             12                  9   -0.122724       -0.137732           0.015008 0.165807 0.580323
                   Basal             12                  9    0.003471       -0.068035           0.071506 0.188598 0.586750
             Macrophages             12                 10   -0.021510       -0.072613           0.051103 0.276605 0.704085
                     AT2             12                 10   -0.052277       -0.063375           0.011098 0.306762 0.715778
```