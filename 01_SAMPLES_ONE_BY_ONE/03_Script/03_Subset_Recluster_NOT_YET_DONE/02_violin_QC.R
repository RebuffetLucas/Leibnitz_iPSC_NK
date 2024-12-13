## @knitr violin_QC

cat(" \n \n")
cat("#### Violin QC ")
cat(" \n \n")

dfCont<-data.frame(
 nFeature_RNA= myObjectSeurat[["nFeature_RNA"]],
 nCount_RNA= myObjectSeurat[["nCount_RNA"]],
 percent.mito= myObjectSeurat[["percent.mito"]],
 orig.ident= myObjectSeurat[["orig.ident"]],
 percent.ribo= myObjectSeurat[["percent.ribo"]]
)

p1<-ggplot(dfCont,aes(x=orig.ident,y=nFeature_RNA,fill=orig.ident))+ 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        .width = 0, 
        justification = -.3, 
        point_colour = NA) + 
    geom_boxplot(
        width = .15, 
        outlier.shape = NA
    ) +
    ## add justified jitter from the {gghalves} package
    gghalves::geom_half_point(
        ## draw jitter on the left
        side = "l", 
        ## control range of jitter
        range_scale = .6, 
        ## add some transparency
        alpha = .1,size=0.1
    )+scale_fill_manual(values = my_palette)

p2<-ggplot(dfCont,aes(x=orig.ident,y=nCount_RNA,fill=orig.ident))+ 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        .width = 0, 
        justification = -.3, 
        point_colour = NA) + 
    geom_boxplot(
        width = .15, 
        outlier.shape = NA
    ) +
    ## add justified jitter from the {gghalves} package
    gghalves::geom_half_point(
        ## draw jitter on the left
        side = "l", 
        ## control range of jitter
        range_scale = .6, 
        ## add some transparency
        alpha = .1,size=0.1
    )+scale_fill_manual(values = my_palette)
p3<-ggplot(dfCont,aes(x=orig.ident,y=percent.mito,fill=orig.ident))+ 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        .width = 0, 
        justification = -.3, 
        point_colour = NA) + 
    geom_boxplot(
        width = .15, 
        outlier.shape = NA
    ) +
    ## add justified jitter from the {gghalves} package
    gghalves::geom_half_point(
        ## draw jitter on the left
        side = "l", 
        ## control range of jitter
        range_scale = .6, 
        ## add some transparency
        alpha = .1,size=0.1
    )+scale_fill_manual(values = my_palette)
p4<-ggplot(dfCont,aes(x=orig.ident,y=percent.ribo,fill=orig.ident))+ 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        .width = 0, 
        justification = -.3, 
        point_colour = NA) + 
    geom_boxplot(
        width = .15, 
        outlier.shape = NA
    ) +
    ## add justified jitter from the {gghalves} package
    gghalves::geom_half_point(
        ## draw jitter on the left
        side = "l", 
        ## control range of jitter
        range_scale = .6, 
        ## add some transparency
        alpha = .1,size=0.1
    )+scale_fill_manual(values = my_palette)


print( p1+p2+p3+p4 )


## @knitr FeatureScatter

cat(" \n \n")
cat("#### Feature Scatter ")
cat(" \n \n")


cat(" \n \n")
print( FeatureScatter(myObjectSeurat, feature1 = "nCount_RNA", feature2 = "percent.mito",group.by = "orig.ident",cols = my_palette,pt.size = 0.5)  + geom_hline(yintercept = VISUAL_TRESHOLD_MITO, linetype = "dashed", color = "red"))
cat(" \n \n")

cat(" \n \n")
print( FeatureScatter(myObjectSeurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident",cols = my_palette,pt.size = 0.5))
cat(" \n \n")

cat(" \n \n")
print( FeatureScatter( myObjectSeurat, feature1 = 'percent.mito', feature2='percent.ribo',group.by = "orig.ident",cols = my_palette,pt.size = 0.5) + geom_vline(xintercept = VISUAL_TRESHOLD_MITO, linetype = "dashed", color = "red"))
cat(" \n \n")



############################################################
####### Same thing with SingleR naming as colour ###########
############################################################

## @knitr FeatureScatter2

Nb_pops=length(unique(myObjectSeurat@meta.data[["BlueprintEncodeData.fine"]]))
mypalette2 = hue_pal()(Nb_pops)
names(mypalette2) <- unique(myObjectSeurat@meta.data[["BlueprintEncodeData.fine"]])
# Change the color of "NK cells" to #FF0000
mypalette2["NK cells"] <- "#FF0000"



cat(" \n \n")
cat("#### Feature Scatter with SingleR Id")
cat(" \n \n")


cat(" \n \n")
print( FeatureScatter(myObjectSeurat, feature1 = "nCount_RNA", feature2 = "percent.mito",group.by = "BlueprintEncodeData.fine",cols = mypalette2, pt.size = 0.5)  + geom_hline(yintercept = VISUAL_TRESHOLD_MITO, linetype = "dashed", color = "red") + geom_vline(xintercept = VISUAL_TRESHOLD_nCOUNT, linetype = "dashed", color = "red"))
cat(" \n \n")


cat(" \n \n")
print( FeatureScatter(myObjectSeurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "BlueprintEncodeData.fine",cols = mypalette2, pt.size = 0.5) + geom_hline(yintercept = VISUAL_TRESHOLD_nFEATURE, linetype = "dashed", color = "red") + geom_vline(xintercept = VISUAL_TRESHOLD_nCOUNT, linetype = "dashed", color = "red"))
cat(" \n \n")


cat(" \n \n")
print( FeatureScatter( myObjectSeurat, feature1 = 'percent.mito', feature2='percent.ribo',group.by = "BlueprintEncodeData.fine",cols = mypalette2, pt.size = 0.5) + geom_vline(xintercept = VISUAL_TRESHOLD_MITO, linetype = "dashed", color = "red"))
cat(" \n \n")


