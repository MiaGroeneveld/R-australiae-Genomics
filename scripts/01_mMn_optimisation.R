# radstackshelpr_visualisation (DeRaad, 2021)
library(RADstackshelpR)

# optimize_m function will generate summary stats on iterative runs
# input
m.out<-optimize_m(m3="m3.vcf",
           m4="m4.vcf",
           m5="m5.vcf",
           m6="m6.vcf",
           m7="m7.vcf")

# Assigning the output of this function to the variable 'm.out' should generate a list containing five objects of class 'data.frame' with the following characteristics:
# 'depth' showing depth per sample for each m value,
# 'snp' showing the number of non-missing SNPs retained in each sample at each m value,
# 'loci' showing the number of non-missing loci retained in each sample at each m value,
# 'snp.R80' showing the total number of SNPs retained at an 80% completeness cutoff,
# 'loci.R80' showing the total number of polymorphic loci retained at an 80% completeness cutoff.

#Use this output list as input for this function, to visualise the effect of varying m on the depth of each sample
vis_depth(output = m.out)

#visualise the effect of varying m on the number of SNPs retained
vis_snps(output = m.out, stacks_param = "m")

#visualise the effect of varying m on the number of loci retained
vis_loci(output = m.out, stacks_param = "m")

#optimize_bigM function will generate summary stats on your 8 iterative runs
bigM.out<-optimize_bigM(M1="bigM1.vcf",
           M2="bigM2.vcf",
           M3="bigM3.vcf",
           M4="bigM4.vcf",
           M5="bigM5.vcf",
           M6="bigM6.vcf",
           M7="bigM7.vcf",
           M8="bigM8.vcf")

#use this function to visualize the effect of varying 'M' on the number of SNPs retained
vis_snps(output = bigM.out, stacks_param = "M")

#visualize the effect of varying 'M' on the number of polymorphic loci retained
vis_loci(output = bigM.out, stacks_param = "M")

#optimize n
n.out<-optimize_n(nequalsMminus1="n2.vcf",
           nequalsM="n3.vcf",
           nequalsMplus1="n4.vcf")

##Assigning the output of this function to the variable 'n.out' should generate a single object of class 'data.frame' showing
#the number of SNPs and loci retained across filtering levels for each value of n.

#visualise the effect of varying n on the number of SNPs retained
vis_snps(output = n.out, stacks_param = "n")

#visualise the effect of varying n on the number of polymorphic loci retained
vis_loci(output = n.out, stacks_param = "n")

#load gridExtra package to combine ggplot visualisations
library(gridExtra)

#combine all of these prior visualisations in a single list
gl<-list()
gl[[1]]<-vis_depth(output = m.out)
# Visualize how different values of m affect average depth in each sample"

gl[[2]]<-vis_snps(output = m.out, stacks_param = "m")
# Visualize how different values of m affect number of SNPs retained.

gl[[3]]<-vis_loci(output = m.out, stacks_param = "m")
# Visualize how different values of m affect number of polymorphic loci retained.

gl[[4]]<-vis_snps(output = bigM.out, stacks_param = "M")
# Visualize how different values of M affect number of SNPs retained.

gl[[5]]<-vis_loci(output = bigM.out, stacks_param = "M")
# Visualize how different values of M affect number of polymorphic loci retained.

gl[[6]]<-vis_snps(output = n.out, stacks_param = "n")
# Visualize how different values of n affect number of SNPs retained.

gl[[7]]<-vis_loci(output = n.out, stacks_param = "n")
# Visualize how different values of n affect number of polymorphic loci retained.

#visualize each item of the list as part of a single grid
grid.arrange(grobs = gl, widths = c(1,1,1,1,1,1),
  layout_matrix = rbind(c(1,1,2,2,3,3),
                        c(4,4,4,5,5,5),
                        c(6,6,6,7,7,7))
)
