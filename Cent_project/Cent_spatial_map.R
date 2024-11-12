# Required packages
# Loading worldclim/cimp6 bioclimatic data
library(terra)
library(geodata)
# displaying images and maps
library(fields)
library(maps)
# Adjusting genotype-environment association models
library(LEA)

# default timeout option is 60s -- increase to 300s
options(timeout = max(300, getOption("timeout")))
# download sample genotypes in the working directory (54.4 MB)
url = "http://membres-timc.imag.fr/Olivier.Francois/Arabidopsis/A_thaliana_chr1.geno"
download.file(url = url, destfile = "./A_thaliana_chr1.geno")
# download sample coordinates in the working directory

url = "http://membres-timc.imag.fr/Olivier.Francois/Arabidopsis/at_coord.coord"
download.file(url = url, destfile = "./at_coord.coord")

genotype = LEA::read.geno("./A_thaliana_chr1.geno")
coordinates = as.matrix(read.table("./at_coord.coord"))

plot(coordinates, cex = .4, col = "darkblue",
     xlab = "Longitude", ylab = "Latitude",
     main = "Sample coordinates", las = 1)
maps::map(add = TRUE, interior = FALSE, col = "grey40")



# Download global bioclimatic data from worldclim
climate <- geodata::worldclim_global(var = 'bio',
                                     res = 10,
                                     download = TRUE,
                                     path=tempdir())

# Download future climate scenario from 'ACCESS-ESM1-5' climate model.
climate_future <- geodata::cmip6_world(model='ACCESS-ESM1-5',
                                       ssp='245',
                                       time='2041-2060',
                                       var='bioc',
                                       download = TRUE,
                                       res=10,
                                       path=tempdir())
# extracting historical environmental data for A. thaliana samples
X.env = terra::extract(x = climate,
                       y = data.frame(coordinates),
                       cells = FALSE)
# remove IDs
X.env = X.env[,-1]
# extracting future environmental data for A. thaliana samples
#X.env_fut = terra::extract(x = climate_future, y = data.frame(coordinates), cells=FALSE)
#X.env_fut = X.env_fut[,-1]


# latent factor GEA model
mod_lfmm = LEA::lfmm2(input = genotype,
                      env = scale(X.env),
                      K = 5,
                      effect.sizes = TRUE)
# get environmental effect sizes
B <- mod_lfmm@B

pv = lfmm2.test(mod_lfmm,
                input = genotype,
                env = scale(X.env),
                full = TRUE)
plot(-log10(pv$pvalue),
     xlab = "SNPs",
     cex = .3, pch = 19, col = "blue")

# define candidate loci for GO analysis
candidates = -log10(pv$pvalue) > 5
# taking all loci for GO analysis
# candidates = -log10(pv$pvalue) > 0
# how many candidate loci?
sum(candidates)



