## ----setup-options, echo = FALSE, results = "hide", message = FALSE------
knitr::opts_chunk$set(comment=NA, fig.align = "center",
                      out.width = "0.7\\linewidth",
                      echo = TRUE, message = TRUE, warning = TRUE,
                      cache = TRUE)
knitr::knit_hooks$set(crop.plot = knitr::hook_pdfcrop)

## ----packages, echo = FALSE, results = "hide", message = FALSE-----------
library("vegan")
data(varespec)
data(varechem)

## ----cca-model-----------------------------------------------------------
cca1 <- cca(varespec ~ ., data = varechem)
cca1

## ----rda-model-----------------------------------------------------------
rda1 <- rda(varespec ~ ., data = varechem)
rda1

## ----eigenvals-----------------------------------------------------------
eigenvals(cca1)

## ----scores--------------------------------------------------------------
str(scores(cca1, choices = 1:4, display = c("species","sites")), max = 1)
head(scores(cca1, choices = 1:2, display = "sites"))

## ----scaling-example, results = "hide"-----------------------------------
scores(cca1, choices = 1:2, display = "species", scaling = 3)

## ----partial-ordination--------------------------------------------------
pcca <- cca(X = varespec,
            Y = varechem[, "Ca", drop = FALSE],
            Z = varechem[, "pH", drop = FALSE])
pcca <- cca(varespec ~ Ca + Condition(pH), data = varechem) ## easier!

## ----triplot-1, fig.height = 5, crop.plot = TRUE, out.width = "0.5\\linewidth"----
plot(cca1)

## ----triplot-2, fig.height = 5, crop.plot = TRUE, out.width = "0.5\\linewidth"----
plot(cca1)

## ----cca-model-build1----------------------------------------------------
vare.cca <- cca(varespec ~ Al + P*(K + Baresoil), data = varechem)
vare.cca

## ----stepwise-1----------------------------------------------------------
upr <- cca(varespec ~ ., data = varechem)
lwr <- cca(varespec ~ 1, data = varechem)
set.seed(1)
mods <- ordistep(lwr, scope = formula(upr), trace = 0)

## ----stepwise-cca--------------------------------------------------------
mods

## ----stepwise-anova------------------------------------------------------
mods$anova

## ----stepwise-reverse----------------------------------------------------
mods2 <- step(upr, scope = list(lower = formula(lwr), upper = formula(upr)), trace = 0,
              test = "perm")
mods2

## ----cca-anova-----------------------------------------------------------
set.seed(42)
(perm <- anova(cca1))

## ----cca-anova-2---------------------------------------------------------
perm

## ----anova-args----------------------------------------------------------
args(anova.cca)

## ----anova-by-axis-------------------------------------------------------
set.seed(1)
anova(mods, by = "axis")

## ----anova-by-term-------------------------------------------------------
set.seed(5)
anova(mods, by = "terms")

## ----anova-by-margin-----------------------------------------------------
set.seed(10)
anova(mods, by = "margin")

## ----meadows-setup-------------------------------------------------------
## load vegan
library("vegan")

## load the data
spp <- read.csv("data/meadow-spp.csv", header = TRUE, row.names = 1)
env <- read.csv("data/meadow-env.csv", header = TRUE, row.names = 1)

## ----meadows-cca-full----------------------------------------------------
m1 <- cca(spp ~ ., data = env)
set.seed(32)
anova(m1)

## ----meadows-cca-full-triplot, fig.height = 5, crop.plot = TRUE, out.width = "0.5\\linewidth"----
plot(m1)

## ----meadows-cca-stepwise------------------------------------------------
set.seed(67)
lwr <- cca(spp ~ 1, data = env)
m2 <- ordistep(lwr, scope = formula(m1), trace = FALSE)
m2

## ----meadows-cca-reduced-triplot, fig.height = 5, crop.plot = TRUE, out.width = "0.5\\linewidth"----
plot(m2)

## ----meadows-cca-anova---------------------------------------------------
m2$anova

## ----meadows-rda---------------------------------------------------------
spph <- decostand(spp, method = "hellinger")
m3 <- rda(spph ~ ., data = env)
lwr <- rda(spph ~ 1, data = env)
m4 <- ordistep(lwr, scope = formula(m3), trace = FALSE)

## ----meadows-rda-reduced-triplot, fig.height = 5, crop.plot = TRUE, out.width = "0.5\\linewidth"----
plot(m4)

## ----meadows-rda-adjrsquare----------------------------------------------
m5 <- ordiR2step(lwr, scope = formula(m3), trace = FALSE)
m5$anova

## ----goodness------------------------------------------------------------
head(goodness(mods))
head(goodness(mods, summarize = TRUE))

## ----inertcomp-----------------------------------------------------------
head(inertcomp(mods, proportional = TRUE))

## ----spenvcor------------------------------------------------------------
spenvcor(mods)

## ----intersetcor---------------------------------------------------------
intersetcor(mods)

## ----shuffle-time-series-------------------------------------------------
shuffle(10, control = how(within = Within(type = "series")))

## ----set-up-toroidal, include = FALSE------------------------------------
set.seed(4)
h <- how(within = Within(type = "grid", ncol = 3, nrow = 3))
perm <- shuffle(9, control = h)

## ----show-toroidal-------------------------------------------------------
matrix(perm, ncol = 3)

## ------------------------------------------------------------------------
plt <- gl(3, 10)
h <- how(within = Within(type = "series"), plots = Plots(strata = plt))

## ----helper-funs---------------------------------------------------------
args(Within)
args(Plots)

## ----how-args------------------------------------------------------------
args(how)

## ----ts-perm-example1----------------------------------------------------
plt <- gl(3, 10)
h <- how(within = Within(type = "series"), plots = Plots(strata = plt))
set.seed(4)
p <- shuffle(30, control = h)
do.call("rbind", split(p, plt)) ## look at perms in context

## ----ts-perm-example2----------------------------------------------------
plt <- gl(3, 10)
h <- how(within = Within(type = "series", constant = TRUE), plots = Plots(strata = plt))
set.seed(4)
p <- shuffle(30, control = h)
do.call("rbind", split(p, plt)) ## look at perms in context

## ----worked-example-devel-1----------------------------------------------
## Analyse the Ohraz data Case study 5 of Leps & Smilauer

## load vegan
library("vegan")

## load the data
spp <- read.csv("data/ohraz-spp.csv", header = TRUE, row.names = 1)
env <- read.csv("data/ohraz-env.csv", header = TRUE, row.names = 1)
molinia <- spp[, 1]
spp <- spp[, -1]

## Year as numeric
env <- transform(env, year = as.numeric(as.character(year)))

## ----worked-example-devel-2----------------------------------------------
## hypothesis 1
c1 <- rda(spp ~ year + year:mowing + year:fertilizer +
          year:removal + Condition(plotid), data = env)

h <- how(within = Within(type = "none"),
         plots = Plots(strata = env$plotid, type = "free"))
set.seed(42)

anova(c1, permutations = h, model = "reduced")
anova(c1, permutations = h, model = "reduced", by = "axis")

## ----worked-example-devel-3----------------------------------------------
## hypothesis 2
c2 <- rda(spp ~ year:mowing + year:fertilizer + year:removal +
          Condition(year + plotid), data = env)
anova(c2, permutations = h, model = "reduced")
anova(c2, permutations = h, model = "reduced", by = "axis")

