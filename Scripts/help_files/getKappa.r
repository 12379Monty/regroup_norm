function (Clustering.v, Labels.v) 
{
    require(concord)
    if (!is.null(names(Clustering.v))) 
        Labels.v <- Labels.v[names(Clustering.v)]
    MS_1.v <- as.numeric(as.character(factor(Labels.v, level = c("MSS", 
        "MSI"), labels = c(1, 2))))
    MS_2.v <- as.numeric(as.character(factor(Labels.v, level = c("MSS", 
        "MSI"), labels = c(2, 1))))
    kappa.1 <- cohen.kappa(cbind(MS_1.v, Clustering.v))
    kappa.2 <- cohen.kappa(cbind(MS_2.v, Clustering.v))
    max(kappa.1$kappa.c, kappa.2$kappa.c)
    kappa.v <- round(max(kappa.1$kappa.c, kappa.2$kappa.c), 2)
}
