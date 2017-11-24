incprops(PrevH = c(0.20), RSE_PrevH = c(0.028),
         PrevR = c(0.10), RSE_PrevR = c(0.094),
         BS_Count = 100000, Boot = TRUE, BMest = "same.test", MDRI = c(200),
         RSE_MDRI = c(0.05), FRR = c(0.01), RSE_FRR = c(0.2),
         BigT = 730)

incprops(PrevH = c(0.20,0.37,0.18), RSE_PrevH = c(0.028,0.1,0.022),
         PrevR = c(0.10,0.18,0.12), RSE_PrevR = c(0.094,0.15,0.05),
         BS_Count = 100000, Boot = TRUE, BMest = 'MDRI.FRR.indep', MDRI = c(200,250,180),
         RSE_MDRI = c(0.05,0.1,0.06), FRR = c(0.01,0.03,0.02), RSE_FRR = c(0.2,0.8,0.1),
         BigT = 730)

