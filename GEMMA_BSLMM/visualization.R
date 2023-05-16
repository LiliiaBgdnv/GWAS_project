#Complex trait

library(data.table)

# Load parameters output
# ==============================================================================
params<-fread("output/gemma_bslmm_output_complex.param.txt",header=T,sep="\t", data.table=F)
# ==============================================================================

# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file="results/complex/top1eff_complex.dsv", quote=F, row.names=F, sep="\t")
write.table(top01, file="results/complex/top0.1eff_complex.dsv", quote=F, row.names=F, sep="\t")
write.table(top001, file="results/complex/top0.01eff_complex.dsv", quote=F, row.names=F, sep="\t")

# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file="results/complex/pip01_complex.dsv", quote=F, row.names=F, sep="\t")
write.table(pip10, file="results/complex/pip10_complex.dsv", quote=F, row.names=F, sep="\t")
write.table(pip25, file="results/complex/pip25_complex.dsv", quote=F, row.names=F, sep="\t")
write.table(pip50, file="results/complex/pip50_complex.dsv", quote=F, row.names=F, sep="\t")
# ==============================================================================
# ==============================================================================
# plot variants PIPs across linkage groups/chromosomes
# ==============================================================================

# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$rs),]

# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))
# ------------------------------------------------------------------------------

# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
png(file="results/complex/pip_plot_complex.png", width=10,height=6,units="in",res=200)

# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0, 0.1),ylab="PIP",xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups
# ------------------------------------------------------------------------------
start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-"light grey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}
# Add variants outside linkage groups
chrs<-c(chrs,"NA")
size<-nrow(params.sort[params.sort$chr=="NA",])
lab.pos<-c(lab.pos, start+size/2)
# ------------------------------------------------------------------------------

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for all variants
# ------------------------------------------------------------------------------
# rank of variants across linkage groups
x<-seq(1,length(params.sort$gamma),1)
# PIP 
y<-params.sort$gamma
# sparse effect size, used for dot size
z<-params.sort$eff
# log-transform to enhance visibility
z[z==0]<-0.00000000001
z<-1/(abs(log(z))*10)
# plot
symbols(x,y,circles=z, bg="black",inches=1/100, fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# highligh high PIP variants (PIP>=0.25)
# ------------------------------------------------------------------------------
# plot threshold line
abline(h=0.05,lty=3,col="dark grey")
# rank of high PIP variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma>=0.01],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma>=0.01]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma>=0.01]
z<-1/abs(log(z))

symbols(x,y,circles=z, bg="red",inches=1/10,fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# add label high PIP variants
text(x,y,labels=params.sort$rs[params.sort$gamma>=0.01], adj=c(0,0), cex=0.8)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# close device
dev.off()
# ==============================================================================
# ==============================================================================
# ==============================================================================


#Simple trait

# Load parameters output
# ==============================================================================
params<-fread("output/gemma_bslmm_output_simple.param.txt",header=T,sep="\t", data.table=F)
# ==============================================================================

# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file="results/simple/top1eff_simple.dsv", quote=F, row.names=F, sep="\t")
write.table(top01, file="results/simple/top0.1eff_simple.dsv", quote=F, row.names=F, sep="\t")
write.table(top001, file="results/simple/top0.01eff_simple.dsv", quote=F, row.names=F, sep="\t")

# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file="results/simple/pip01_simple.dsv", quote=F, row.names=F, sep="\t")
write.table(pip10, file="results/simple/pip10_simple.dsv", quote=F, row.names=F, sep="\t")
write.table(pip25, file="results/simple/pip25_simple.dsv", quote=F, row.names=F, sep="\t")
write.table(pip50, file="results/simple/pip50_simple.dsv", quote=F, row.names=F, sep="\t")
# ==============================================================================
# ==============================================================================
# plot variants PIPs across linkage groups/chromosomes
# ==============================================================================

# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$rs),]

# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))
# ------------------------------------------------------------------------------

# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
png(file="results/simple/pip_plot_simple.png", width=10,height=6,units="in",res=200)

# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0, 0.1),ylab="PIP",xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups
# ------------------------------------------------------------------------------
start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-"light grey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}
# Add variants outside linkage groups
chrs<-c(chrs,"NA")
size<-nrow(params.sort[params.sort$chr=="NA",])
lab.pos<-c(lab.pos, start+size/2)
# ------------------------------------------------------------------------------

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for all variants
# ------------------------------------------------------------------------------
# rank of variants across linkage groups
x<-seq(1,length(params.sort$gamma),1)
# PIP 
y<-params.sort$gamma
# sparse effect size, used for dot size
z<-params.sort$eff
# log-transform to enhance visibility
z[z==0]<-0.00000000001
z<-1/(abs(log(z))*10)
# plot
symbols(x,y,circles=z, bg="black",inches=1/100, fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# highligh high PIP variants (PIP>=0.25)
# ------------------------------------------------------------------------------
# plot threshold line
abline(h=0.05,lty=3,col="dark grey")
# rank of high PIP variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma>=0.01],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma>=0.01]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma>=0.01]
z<-1/abs(log(z))

symbols(x,y,circles=z, bg="red",inches=1/10,fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# add label high PIP variants
text(x,y,labels=params.sort$rs[params.sort$gamma>=0.01], adj=c(0,0), cex=0.8)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# close device
dev.off()
# ==============================================================================
# ==============================================================================
# ==============================================================================

#Complex trait generated data

# Load parameters output
# ==============================================================================
params<-fread("output/gemma_bslmm_output_complex_gen.param.txt",header=T,sep="\t", data.table=F)
# ==============================================================================

# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file="results/complex_gen/top1eff_complex_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(top01, file="results/complex_gen/top0.1eff_complex_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(top001, file="results/complex_gen/top0.01eff_complex_gen.dsv", quote=F, row.names=F, sep="\t")

# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file="results/complex_gen/pip01_complex_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(pip10, file="results/complex_gen/pip10_complex_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(pip25, file="results/complex_gen/pip25_complex_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(pip50, file="results/complex_gen/pip50_complex_gen.dsv", quote=F, row.names=F, sep="\t")
# ==============================================================================
# ==============================================================================
# plot variants PIPs across linkage groups/chromosomes
# ==============================================================================

# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$rs),]

# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))
# ------------------------------------------------------------------------------

# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
png(file="results/complex_gen/pip_plot_complex_gen.png", width=10,height=6,units="in",res=200)

# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0, 0.1),ylab="PIP",xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups
# ------------------------------------------------------------------------------
start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-"light grey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}
# Add variants outside linkage groups
chrs<-c(chrs,"NA")
size<-nrow(params.sort[params.sort$chr=="NA",])
lab.pos<-c(lab.pos, start+size/2)
# ------------------------------------------------------------------------------

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for all variants
# ------------------------------------------------------------------------------
# rank of variants across linkage groups
x<-seq(1,length(params.sort$gamma),1)
# PIP 
y<-params.sort$gamma
# sparse effect size, used for dot size
z<-params.sort$eff
# log-transform to enhance visibility
z[z==0]<-0.00000000001
z<-1/(abs(log(z))*10)
# plot
symbols(x,y,circles=z, bg="black",inches=1/100, fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# highligh high PIP variants (PIP>=0.25)
# ------------------------------------------------------------------------------
# plot threshold line
abline(h=0.05,lty=3,col="dark grey")
# rank of high PIP variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma>=0.01],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma>=0.01]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma>=0.01]
z<-1/abs(log(z))

symbols(x,y,circles=z, bg="red",inches=1/10,fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# add label high PIP variants
text(x,y,labels=params.sort$rs[params.sort$gamma>=0.01], adj=c(0,0), cex=0.8)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# close device
dev.off()
# ==============================================================================
# ==============================================================================
# ==============================================================================

#Simple trait generated data

# Load parameters output
# ==============================================================================
params<-fread("output/gemma_bslmm_output_simple_gen.param.txt",header=T,sep="\t", data.table=F)
# ==============================================================================

# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file="results/simple_gen/top1eff_simple_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(top01, file="results/simple_gen/top0.1eff_simple_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(top001, file="results/simple_gen/top0.01eff_simple_gen.dsv", quote=F, row.names=F, sep="\t")

# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file="results/simple_gen/pip01_simple_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(pip10, file="results/simple_gen/pip10_simple_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(pip25, file="results/simple_gen/pip25_simple_gen.dsv", quote=F, row.names=F, sep="\t")
write.table(pip50, file="results/simple_gen/pip50_simple_gen.dsv", quote=F, row.names=F, sep="\t")
# ==============================================================================
# ==============================================================================
# plot variants PIPs across linkage groups/chromosomes
# ==============================================================================

# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$rs),]

# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))
# ------------------------------------------------------------------------------

# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
png(file="results/simple_gen/pip_plot_simple_gen.png", width=10,height=6,units="in",res=200)

# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0, 0.1),ylab="PIP",xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups
# ------------------------------------------------------------------------------
start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-"light grey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}
# Add variants outside linkage groups
chrs<-c(chrs,"NA")
size<-nrow(params.sort[params.sort$chr=="NA",])
lab.pos<-c(lab.pos, start+size/2)
# ------------------------------------------------------------------------------

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for all variants
# ------------------------------------------------------------------------------
# rank of variants across linkage groups
x<-seq(1,length(params.sort$gamma),1)
# PIP 
y<-params.sort$gamma
# sparse effect size, used for dot size
z<-params.sort$eff
# log-transform to enhance visibility
z[z==0]<-0.00000000001
z<-1/(abs(log(z))*10)
# plot
symbols(x,y,circles=z, bg="black",inches=1/100, fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# highligh high PIP variants (PIP>=0.25)
# ------------------------------------------------------------------------------
# plot threshold line
abline(h=0.05,lty=3,col="dark grey")
# rank of high PIP variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma>=0.01],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma>=0.01]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma>=0.01]
z<-1/abs(log(z))

symbols(x,y,circles=z, bg="red",inches=1/10,fg=NULL,add=T, lables = NULL)
# ------------------------------------------------------------------------------

# add label high PIP variants
text(x,y,labels=params.sort$rs[params.sort$gamma>=0.01], adj=c(0,0), cex=0.8)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# close device
dev.off()
# ==============================================================================
