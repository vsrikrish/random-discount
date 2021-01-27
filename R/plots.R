library(reshape2)
library(ggplot2)
source('R/random_discount.R')
nsamp <- 1e5
nyr <- 400

# load data
dat <- readRDS('data/discount.rds')
# remove NAs
dat <- dat[!is.na(dat$DR),]
dat$Year <- as.numeric(dat$Year)

models <- c('random', 'mean', 'drift')
parnames <- vector('list', 3)
names(parnames) <- models
parnames[['random']] <- c('rho1', 'rho2', 'rho3', 'sigma_sq')
parnames[['mean']] <- c('eta', 'rho1', 'rho2', 'rho3', 'sigma_sq')
parnames[['drift']] <- c('a', 'b', 'rho1', 'rho2', 'rho3', 'sigma_sq')

m <- lapply(models, function(n) readRDS(paste0('output/mcmc-', n, '.rds')))
names(m) <- models
pars <- lapply(m, function(l) l$samples[sample(1:nrow(l$samples), nsamp, replace=TRUE),])
rt <- mapply(function(p, n, pn) {match.fun(paste0(n, '_discount'))(log(dat$DR), pars=p, N=nyr, parnames=pn)}, pars, models, parnames, USE.NAMES=TRUE, SIMPLIFY=FALSE)
df <- lapply(rt, compute_factor)

ce <- lapply(df, ce_discount)
condf <- lapply(c(4/100, 1.5/100), compute_factor, t=nyr)
names(condf) <- c('const4', 'const1.5')
ce <- c(ce, condf)
cem <- melt(ce)
cem$year <- rep(0:nyr, length(ce))

# plot certainty-equivalent rates to get value of $100 today
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
keylabels <- c('Constant 1.5%/yr (OMB)', 'Constant 4%/yr', 'Drift', 'Mean-Reverting', 'Random Walk')
p <- ggplot(cem) + geom_line(aes(x=year, y=value*100, color=L1, linetype=L1), size=1.25) + 
  scale_color_manual('Discount Model', values=cbbPalette[c(1, 1, 6, 4, 7)], labels=keylabels) +
  scale_linetype_manual('Discount Model', values=c(2, 4, 1, 1, 1), labels=keylabels) +
  scale_y_log10('Certainty-Equivalent \nValue Today of $100', labels=scales::dollar_format(), expand=c(0, 0), limits=c(1e-2, 100)) +
  scale_x_continuous('Years From 2019', expand=c(0.01, 0.01)) +
  theme_classic(base_size=16) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position=c(0.01, 0.01), legend.justification=c(0, 0)) +
  annotation_logticks(sides='l') +
png('figures/ce_values.png', res=600, height=6, width=8, units='in')
p
dev.off()

# plot quantiles for discount rates
rt_q <- lapply(rt, function(l) t(apply(l, 2, quantile, probs=c(0.05, 0.5, 0.95))))
condrt <- lapply(c(4, 1.5), function(x) data.frame('5%'=rep(NA, nyr), '50%'=rep(x, nyr), '95%'=rep(NA, nyr), check.names=FALSE))
names(condrt) <- c('const4', 'const1.5')
rt_m <- do.call(rbind, c(condrt, rt_q))
colnames(rt_m) <- c('low', 'median', 'upper')
rt_m$model <- rep(c(names(condrt), names(rt_q)), each=nyr)
rt_m$year <- rep(1:nyr, length(names(condrt)) + length(names(rt_q)))
p <- ggplot(rt_m) + geom_ribbon(aes(x=year, ymin=low/100, ymax=upper/100, fill=model), alpha=0.2) +
  geom_line(aes(x=year, y=median/100, color=model, linetype=model)) +
  geom_line(data=dat, aes(x=Year-2019, y=DR/100), color='grey') +
  scale_color_manual('Discount Model', values=cbbPalette[c(1, 1, 6, 4, 7)], labels=keylabels) +
  scale_fill_manual('Discount Model', values=cbind(cbbPalette[c(1, 1, 6, 4, 7)]), labels=keylabels) +
  scale_linetype_manual('Discount Model', values=c(2, 4, 1, 1, 1), labels=keylabels) +
  scale_x_continuous('Years from 2019', expand=c(0.01, 0.01), limits=c(-100, 400)) +
  scale_y_continuous('Discount Rate (%/yr)', expand=c(0,0), labels=scales::percent_format(), breaks=seq(0, .25, .05)) +
  theme_classic(base_size=16) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position=c(0.01, 0.99), legend.justification=c(0, 1)) +
  guides(fill = guide_legend(override.aes = list(fill=c(NA, NA, cbbPalette[c(6, 4, 7)]))))
png('figures/rates.png', res=600, height=6, width=8, units='in')
p
dev.off()

                    
