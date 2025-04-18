# Data load ---------------------------------------------------------------
library(tidyverse)

theme_set(theme_bw() + theme(
    text = element_text(family = "serif", size = 15),
    legend.position = "bottom")
)
manual.aic <- function(a){
    if(!(a$family$family %in% c("quasipoisson", "poisson"))){return(NA)}
    loglik <- sum(dpois(a$model[,1], fitted.values(a), log = TRUE))
    phi <- summary(a)$dispersion
    -2 * loglik + 2 * summary(a)$df[3] * phi
}

corr <- function(df, x, y, grp){
    # x = "abu"; # y = "shan"; # grp = c("year", "zone")
    # unify arguments 
    if(is.numeric(x)){x <- colnames(df)[x]}
    if(is.numeric(y)){x <- colnames(df)[y]}
    if(is.numeric(grp)){grp <- colnames(df)[grp]}
    # run
    df %>% 
        unite(all_of(grp), col = type, sep = "_") %>% 
        select(type, all_of(c(x, y))) %>% 
        rename(x = 2, y = 3) %>% 
        split(.$type) %>% 
        map(~cor.test(.x$x, .x$y, method = "spearman")) %>% 
        map_dfr(~tibble(r = .x$estimate, p = .x$p.value), .id = "type") %>% 
        separate(type, sep = "_", into = grp)
}



long <- readxl::read_excel("Data/Carabidae_25.01.2023.xlsx", 
                           sheet = "main_data", range = "A1:H887") %>% 
    arrange(desc(site)) %>% 
    mutate(tur = factor(tur), 
           site = fct_inorder(site),
           zone = case_when(zone == "fon" ~ "фоновая", 
                            zone == "bufer" ~ "буферная",
                            zone == "superimpact" ~ "суперимпактная", 
                            TRUE ~ "импактная"),
           zone = factor(zone, levels = c("фоновая", "буферная", "импактная", "суперимпактная"))) %>% 
    mutate(km = str_extract(site, "[:digit:]{1,}"),
           km = as.numeric(km), 
           year = as.factor(year), 
           .after = "site")
wide <- long %>% 
    select(-num) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    mutate(id = paste0(year, tur, site, plot), .before = year) %>% 
    select(-no_insects)

div <- tibble(wide[,2:7], 
              abu = apply(wide[,8:ncol(wide)], 1, sum),
              nsp = apply(wide[,8:ncol(wide)], 1, function(a){length(a[a>0])}),
              shan= vegan::diversity(wide[,8:ncol(wide)], 
                                     MARGIN = 1, index = "shannon", base = exp(1))
              # dom = apply(wide[,7:ncol(wide)], 1, function(a){a <- a/sum(a); return(max(a))})
)

div2 <- long %>% # div2 - turs are united
    filter(taxa != "no_insects") %>% 
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0)
div2 <- tibble(div2[,1:5], 
               abu = apply(div2[,6:ncol(div2)], 1, sum),
               nsp = apply(div2[,6:ncol(div2)], 1, function(a){length(a[a>0])}),
               shan= vegan::diversity(div2[,6:ncol(div2)], 
                                      MARGIN = 1, index = "shannon", base = exp(1))
)
div
# Базовые подсчёты --------------------------------------------------------
cat("total number of species:", length(unique(long$taxa)))

cat("number of species by years:")
long %>% 
    select(taxa,year) %>% 
    distinct() %>% 
    group_by(year) %>% 
    summarise(nsp_total = n(), .groups = "drop") %>% 
    pivot_wider(names_from = year, values_from = nsp_total)

cat("Особи по зонам и годам:")
long %>% 
    select(taxa, year, num) %>% 
    group_by(year) %>% 
    summarise(n_individuals = sum(num), .groups = "drop") %>% 
    pivot_wider(names_from = year, values_from = n_individuals)

cat("Обилие по годам и турам учётов!!!:")
long %>% 
    group_by(year, tur) %>% 
    summarise(abu = mean(abu), .groups = "drop") 

# dominant complex
dom_comp <- long %>% 
    mutate(zone = fct_collapse(zone, `импактная` = c("импактная", "суперимпактная"))) %>% 
    group_by(zone, year, taxa) %>% 
    summarise(abu = sum(abu), .groups = "drop_last") %>% 
    mutate(p = round(abu/sum(abu)*100, 1)) %>% 
    ungroup()
dom_comp <- dom_comp %>% 
    filter(taxa %in% unique(dom_comp$taxa[dom_comp$p >= 5])) %>% 
    # filter(p >= 5) %>% 
    transmute(yz = paste0(year, "_", zone), taxa, p) %>% 
    pivot_wider(names_from = yz, values_from = p, values_fill = 0) %>%
    arrange(taxa)

library(formattable)
formattable(dom_comp, 
            align = "l",
            list(`2009_фоновая` =  color_bar("lightgreen"),
                 `2014_фоновая` =  color_bar("lightgreen"),
                 `2009_буферная` = color_bar("gold"),
                 `2014_буферная` = color_bar("gold"), 
                 `2009_импактная` =color_bar("coral"),
                 `2014_импактная` =color_bar("coral"), 
                 taxa  = formatter("span", style = ~ style(color = "grey60",font.style = "italic")))) 

# Diversity visuzlisation ------------------------------------------------------------
div %>% 
    mutate(type = year #paste0(tur, " тур, ", year), 
           # km = as.factor(km)
           ) %>% 
    # rbind(., mutate(., type = substr(type, nchar(type)-3, nchar(type)))) %>% 
    ggplot(aes(x = km, y = abu,  color = tur)) + 
    geom_point() + 
    scale_color_manual(values = c("cyan", "coral")) +
    facet_wrap(~type, scales = "free", ncol = 2) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Обилие в разные годы и туры") + 
    scale_fill_manual(values = c("white", "grey85", "grey50", "black")) +
    guides(fill="none") 
ggsave("Plot2.svg", width = 16, height = 11, units = "cm")

div %>% ##
    mutate(type = year #paste0(tur, " тур, ", year), 
           # km = as.factor(km)
           ) %>% 
    # 
    # rbind(., mutate(., type = substr(type, nchar(type)-3, nchar(type)))) %>% 
    # turn off to get journal picture
    ggplot(aes(x = km, y = nsp, color = tur)) + 
    geom_point() + 
    scale_color_manual(values = c("cyan", "coral")) +
    facet_wrap(~year, scales = "free", ncol = 2) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Количество видов", 
         subtitle = "Видовое богатство в разные годы и туры") + 
    scale_fill_manual(values = c("white", "grey85", "grey50", "black")) +
    guides(fill="none")
ggsave("Plot3.svg", width = 16, height = 11*0.6, units = "cm")

div %>% 
    mutate(type = year #paste0(tur, " тур, ", year),
           # km = as.factor(km)
           ) %>%
    # rbind(., mutate(., type = substr(type, nchar(type)-3, nchar(type)))) %>%
    ggplot(aes(x = km, y = shan, color = tur)) + 
    geom_point() + 
    scale_color_manual(values = c("cyan", "coral")) +
    facet_wrap(~year, scales = "free", ncol = 2) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Индекс Шеннона", 
         subtitle = "Таксономическое разнообразие в разные годы и туры") + 
    scale_fill_manual(values = c("white", "grey85", "grey50", "black")) +
    guides(fill="none")
ggsave("Plot4.svg", width = 16, height = 11*0.6, units = "cm")

# Abundance models selection ---------------------------------------------------
fits.abu <- expand_grid(
    # yy = "abu",
    ff = c("year", "year + zone", "year + km", "year + km_log"),
    tt = "", #c("tur_united", "", " + tur"), 
    # учитывать туры лучше, чем не учитывать
    # но объединять ещё лучше, чем учитывать
    super_imp = c(""), #"s.i._excl.", 
    # суперимпакт лучше не исключать, с ним r2 лучше (0.4 -> 0.5)
    distr = c("gaussian") # "poisson", "quasipoisson", 
    # if gaussian then log = true
    # poisson is removed: it produces the same r2 with a little bigger AIC
    # but its almost always doesn't pass via overdispersion test
    # quasipoisson is removed: AIC is incredibly high (in comp. with gauss.)
    # and shapiro test for residuals fails
) %>% 
    mutate(
        ff = case_when(tt == " + tur" ~ paste0(ff, tt), TRUE ~ ff), 
        ff = paste0("abu ~ ", ff), 
        tt = case_when(str_detect(tt, "uni") ~ tt, TRUE ~ "")
        ) %>% 
    mutate(nm = paste0(ff, "_", tt, "_", super_imp, "_", distr))

fits.abu <- fits.abu %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){ 
        if(a$tt == "tur_united"){df <- div2} else {df <- div}
        if(a$super_imp == "s.i._excl."){ 
            df <- filter(df, zone != "суперимпактная")
        }
        if(a$distr == "gaussian"){df <- mutate(df, abu = log10(abu+1))}
        ff <- a$ff
        if(str_detect(ff, "km_log")){
            ff <- str_replace_all(ff, "km_log", "km")
            df <- mutate(df, km = log10(km))
        }
        if(str_detect(ff, "km_sqrt")){
            ff <- str_replace_all(ff, "km_sqrt", "km")
            df <- mutate(df, km = sqrt(km))
        }
        glm(formula = ff, data = df, family = a$distr)
    }) %>% 
    `names<-`(fits.abu$nm) %>% 
    mutate(fits.abu, fit = ., info = lapply(., summary))

fits.abu <- fits.abu %>% 
    pull(fit) %>% 
    lapply(function(a){ 
        resid = residuals(a)
        ftted = fitted.values(a)
        if(a$family$family == "poisson"){
            suppressMessages(overdispersion <- a |> 
                                 performance::check_overdispersion() |>
                                 capture.output() |> 
                                 str_subset("value") |> 
                                 str_replace_all("[:alpha:]|[:blank:]|[:symbol:]|-", "") |> 
                                 as.numeric())
        } else {overdispersion <- NA}
        
        tibble(aic = AIC(a), 
               aic2 = manual.aic(a),
               shapiro = shapiro.test(resid)$p.value,
               overdispersion = overdispersion,
               var.res = var(resid), 
               var.ftd = var(ftted)
        ) %>% 
            mutate(r2 = var.ftd/(var.res + var.ftd)) %>% 
            return()
    }) %>% 
    map_df(rbind, .id = "nm") %>% 
    left_join(fits.abu, ., by = "nm") %>% 
    mutate(aic = case_when(is.na(aic) ~ aic2, TRUE ~ aic), .keep = "unused")

list(
    lm = lm(abu ~ km, data = mutate(div2, abu = log(abu+1))), 
    x2 = lm(abu ~ km + km2, data = mutate(div2, abu = log(abu+1), km2 = km^2)),
    x12= lm(abu ~ km + km12, data = mutate(div2, abu = log(abu+1), km12 = sqrt(km))),
    seg= segmented::segmented(lm(abu ~ km, data = mutate(div2, abu = log(abu+1))), seg.Z = ~km)
) %>% 
    map(~tibble(aic = AIC(.x), r2 = summary(.x)$r.squared))
    
summary(lm_shan)

summary(x2_shan)
seg_shan <- 
summary(seg_shan)

fits.abu %>% select(-fit, -info, -var.res, -var.ftd, -overdispersion) %>% writexl::write_xlsx("tmp.xlsx")

fits.abu %>% 
    filter(str_detect(ff, "site")) %>%
    filter(is.na(overdispersion) | overdispersion > 0.05, 
           str_detect(ff, "tur"), 
           distr == "quasipoisson") %>% 
    View

# Abundance model ---------------------------------------------------------
div %>% group_by(year, tur) %>% summarise(a = mean(abu))

abu_0 <- glm(abu~ year * tur + site, 
             family = "quasipoisson",
             data = filter(div,  zone != "суперимпактная"))
abu_si<- glm(abu~ year * tur + site, 
             family = "quasipoisson", 
             data = div)
summary(abu_0)
summary(abu_si)

cat("AIC_0 =",  manual.aic(abu_0))
cat("AIC_si =", manual.aic(abu_si))
with(summary(abu_0), 1 - deviance/null.deviance) |> round(3) |> cat("model_0: r^2 =", ... = _)
with(summary(abu_si), 1 - deviance/null.deviance)|> round(3) |> cat("model_si: r^2 =", ... = _)

TukeyHSD(aov(abu_0))
TukeyHSD(aov(abu_si))

list(a1 = summary(abu_0)$coefficients, 
     a2 = data.frame(r2 = with(summary(abu_0), 1 - deviance/null.deviance)), 
     a3 = data.frame(aic = manual.aic(abu_0)),
     b1 = summary(abu_si)$coefficients, 
     b2 = data.frame(r2 = with(summary(abu_si), 1 - deviance/null.deviance)), 
     b3 = data.frame(aic = manual.aic(abu_si))
) %>% 
    lapply(as.data.frame) %>% 
    lapply(rownames_to_column) %>% 
    writexl::write_xlsx("abundance.result.xlsx")

# Для оценки значимости категориального предиктора в целом, 
# а не по отдельным градациям, вы можете использовать тест отношения правдоподобия. 
# Для этого вам необходимо подогнать две модели: 
#       одну с категориальным предиктором и другую без него. 
# Затем вы можете использовать функцию anova() в R для сравнения этих двух моделей 
# и проверки значимости категориального предиктора. Например:
#       
#       model1 <- glm(y ~ x1 + x2 + factor, data = mydata)
#       model2 <- glm(y ~ x1 + x2, data = mydata)
#       anova(model1, model2, test = "LRT")
# В этом примере y - это зависимая переменная, 
# x1 и x2 - количественные предикторы, а
# factor - категориальный предиктор с 5 градациями фактора. 
# Функция anova() сравнивает две модели 
# и выполняет тест отношения правдоподобия (LRT) 
# для проверки значимости категориального предиктора1.

# correlation -------------------------------------------------------------
corr.data <- wide %>% 
    select(-id, -tur, -zone, -plot, -site) %>% 
    group_by(year, km) %>% 
    summarise_all(mean) %>% 
    ungroup() %>% 
    mutate(
        nsp = apply(.[3:ncol(.)], 1, function(a){length(a[a>0])}),
        shn = apply(.[3:ncol(.)], 1, function(a){vegan::diversity(a)}),
        .before = 3
    ) %>% 
    select(1:4)

corr.result <- list()

corr.result$abu_10 <- div %>% 
    group_by(year, tur, km) %>% 
    summarise(abu = mean(abu), .groups = "drop") %>% 
    corr("km", "abu", c("year", "tur"))
corr.result$abu_9 <- div %>% 
    filter(km > 1) %>%
    group_by(year, tur, km) %>% 
    summarise(abu = mean(abu), .groups = "drop") %>% 
    corr("km", "abu", c("year", "tur"))
corr.result$nsp_10 <- corr.data %>% 
    corr("km", "nsp", "year")
corr.result$nsp_9 <- corr.data %>% 
    filter(km > 1) %>% 
    corr("km", "nsp", "year")
corr.result$shn_10 <- corr.data %>% 
    corr("km", "shn", "year")
corr.result$shn_9 <- corr.data %>% 
    filter(km > 1) %>% 
    corr("km", "shn", "year")

writexl::write_xlsx(corr.result, "correlation_04.07.2023.xlsx")





cor.sites10 <- div %>% 
    # filter(km > 1) %>% 
    group_by(year, tur, km) %>% 
    summarise(abu = mean(abu), .groups = "drop") %>% 
    unite(type, year, tur) %>% 
    split(.$type) %>% 
    lapply(function(a){
        r <- cor.test(a$km, a$abu, method = "spearman")
        data.frame(est = r$estimate, pval = r$p.value) |> 
            remove_rownames() %>% 
            mutate_all(function(b)round(b, 3)) %>% 
            transmute(lab = paste0("ρ = ", est, "\np.value = ", pval))
    }) %>% 
    map_dfr(rbind, .id = "type")

div %>% 
    group_by(year, tur, km) %>% 
    summarise(abu = mean(abu), 
              .groups = "drop_last") %>% 
    mutate(km2 = rank(km)) %>% 
    ungroup() %>% 
    unite(type, year, tur) %>% 
    ggplot(aes(x = km2, y = abu, color = type)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = "y~x", se = FALSE) + 
    geom_label(aes(x = 3.1, y = 400, label = lab), data = cor.sites10, 
               alpha = 0.6, size = 3) +
    facet_wrap(~type) +
    scale_x_continuous(breaks = 1:10,
                       labels = c(1, 4, 5, 9, 11, 12, 18, 26, 27, 32)) +
    labs(x = "Distance from smelter (not to scale)", 
         y = "Abundance", 
         subtitle = "all plots are united") +
    guides(color=guide_legend(title="years & turs")) +
    theme(panel.grid.minor = element_blank())
ggsave("abu.cor10.png", width = 297*0.7, height = 210*0.7, units = "mm")

cor.sites9 <- div %>% 
    filter(km > 1) %>%
    group_by(year, tur, km) %>% 
    summarise(abu = mean(abu), .groups = "drop") %>% 
    unite(type, year, tur) %>% 
    split(.$type) %>% 
    lapply(function(a){
        r <- cor.test(a$km, a$abu, method = "spearman")
        data.frame(est = r$estimate, pval = r$p.value) |> 
            remove_rownames() %>% 
            mutate_all(function(b)round(b, 3)) %>% 
            transmute(lab = paste0("ρ = ", est, "\np.value = ", pval))
    }) %>% 
    map_dfr(rbind, .id = "type")

div %>% 
    filter(km > 1) %>%
    group_by(year, tur, km) %>% 
    summarise(abu = mean(abu), 
              .groups = "drop_last") %>% 
    mutate(km2 = rank(km)+1) %>% 
    ungroup() %>% 
    unite(type, year, tur) %>% 
    ggplot(aes(x = km2, y = abu, color = type)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = "y~x", se = FALSE) + 
    geom_label(aes(x = 3.8, y = 400, label = lab), data = cor.sites9, 
               alpha = 0.6, size = 3) +
    facet_wrap(~type) +
    scale_x_continuous(breaks = 1:10,
                       labels = c(1, 4, 5, 9, 11, 12, 18, 26, 27, 32)) +
    labs(x = "Distance from smelter (not to scale)", 
         y = "Abundance", 
         subtitle = "all plots are united, superimpact removed") +
    guides(color=guide_legend(title="years & turs")) +
    theme(panel.grid.minor = element_blank())
ggsave("abu.cor10.png", width = 297*0.7, height = 210*0.7, units = "mm")

rbind(cor.sites10, cor.sites9) %>% 
    separate(lab, sep = "\n", into = c("r", "p")) %>% 
    mutate(year = type, type = c(rep("10_sites", 4), rep("9_sites", 4)))

cor1 <- wide %>% 
    select(-id, -tur, -zone, -plot, -site) %>% 
    group_by(year, km) %>% 
    summarise_all(mean) %>% 
    ungroup() %>% 
    mutate(
        nsp = apply(.[3:ncol(.)], 1, function(a){length(a[a>0])}),
        shn = apply(.[3:ncol(.)], 1, function(a){vegan::diversity(a)}),
        .before = 3
    ) %>% 
    select(1:4)


cor1 %>% 
    split(.$year) %>% 
    lapply(function(a){
        r <- cor.test(a$km, a$nsp, method = "spearman")
        data.frame(r = r$estimate, p = r$p.value) |> 
            remove_rownames()
    }) %>% 
    map_dfr(rbind, .id = "year")


# Diversity (number of species) models --------------------------------------------------------
fits.nsp <- expand_grid(
    yy = "nsp",
    xx = c("year + site", "year + tur + site", "year * tur + site"), #c("zone", "site", "km")),
    super_imp = c("S.I.", ""),
    distr = c("poisson", "quasipoisson")
) %>% 
    unite(ff, yy, xx, sep = "~ ") %>%  
    mutate(nm = paste0(ff, "_", super_imp, "_", distr))

fits.nsp <- fits.nsp %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){ 
        if(str_detect(a$ff, "tur")){
            df <- div
        }else{
            df <- div2
        }
        if(a$super_imp != "S.I."){
            df <- filter(df, zone != "суперимпактная")
        }
        glm(formula = a$ff, data = df, family = a$distr)
    }) %>% 
    `names<-`(fits.nsp$nm) %>% 
    mutate(fits.nsp, fit = ., info = lapply(., summary)) 

fits.nsp <- fits.nsp %>% 
    pull(fit) %>% 
    lapply(function(a){ 
        resid = residuals(a)
        ftted = fitted.values(a)
        if(a$family$family == "poisson"){
            suppressMessages(overdispersion <- a |> 
                                 performance::check_overdispersion() |>
                                 capture.output() |> 
                                 str_subset("value") |> 
                                 str_replace_all("[:alpha:]|[:blank:]|[:symbol:]|-", "") |> 
                                 as.numeric())
        } else {overdispersion <- NA}
        
        tibble(aic = AIC(a), 
               aic2 = manual.aic(a),
               shapiro = shapiro.test(resid)$p.value,
               overdispersion = overdispersion,
               var.res = var(resid), 
               var.ftd = var(ftted)
        ) %>% 
            mutate(r2 = var.ftd/(var.res + var.ftd)) %>% 
            return()
    }) %>% 
    map_df(rbind, .id = "nm") %>% 
    left_join(fits.nsp, ., by = "nm") %>% 
    mutate(aic = case_when(is.na(aic) ~ aic2, TRUE ~ aic), .keep = "unused")

fits.nsp <- fits.nsp %>% 
    filter(str_detect(distr, "quasi", negate = TRUE), 
           str_detect(ff, "tur", negate = TRUE)) %>% 
    arrange(super_imp)

list(a1 = fits.nsp$info[[1]]$coefficients, 
     a2 = data.frame(r2 = with(fits.nsp$info[[1]], 1 - deviance/null.deviance)), 
     a3 = data.frame(aic = manual.aic(fits.nsp$fit[[1]])),
     b1 = fits.nsp$info[[2]]$coefficients, 
     b2 = data.frame(r2 = with(fits.nsp$info[[2]], 1 - deviance/null.deviance)), 
     b3 = data.frame(aic = manual.aic(fits.nsp$fit[[2]]))
) %>% 
    lapply(as.data.frame) %>% 
    lapply(rownames_to_column) %>% 
    writexl::write_xlsx("nsp.result.xlsx")

# Diversity (shannon index) models --------------------------------------------------------
fits.shn <- expand_grid(
    yy = "shan",
    xx = "year", #c("zone", "site", "km")),
    zz = c(" + site", " + tur + site", " * tur + site"),
    super_imp = c("S.I.", ""),
    distr =  "gaussian", 
    tursum = c(TRUE, FALSE)) %>% 
    mutate(xx = paste0(xx, zz), .keep = "unused") %>% 
    unite(ff, yy, xx, sep = "~ ") %>% 
    filter(tursum == FALSE | str_detect(ff, "tur", negate = TRUE)) %>% 
    mutate(nm = paste0(ff, "_", super_imp, "_", distr, "..", tursum))

fits.shn <- fits.shn %>% 
    # slice(3) %>%
    split(1:nrow(.)) %>% 
    lapply(function(a){
        if(a$tursum){df <- div2} else {df <- div}
        if(a$super_imp != "S.I."){ 
            df <- filter(df, zone != "суперимпактная")
        }
        glm(formula = a$ff, data = df, family = a$distr)
    }) %>% 
    `names<-`(fits.shn$nm) %>% 
    mutate(fits.shn, fit = ., info = lapply(., summary)) %>% 
    arrange(nm)

fits.shn <- fits.shn %>% 
    pull(fit) %>% 
    lapply(function(a){ 
        resid = residuals(a)
        ftted = fitted.values(a)
        tibble(aic = AIC(a), 
               shapiro = shapiro.test(resid)$p.value,
               var.res = var(resid), 
               var.ftd = var(ftted)
        ) %>% 
            mutate(r2 = var.ftd/(var.res + var.ftd)) %>% 
            return()
    }) %>% 
    map_df(rbind, .id = "nm") %>% 
    left_join(fits.shn, ., by = "nm") %>% 
    arrange(nm)

fits.shn
fits.shn <- fits.shn %>% 
    filter(str_detect(nm, "tur|FALSE", negate = TRUE)) %>% 
    arrange(aic)

list(a1 = fits.shn$info[[1]]$coefficients, 
     a2 = data.frame(r2 = with(fits.shn$info[[1]], 1 - deviance/null.deviance)), 
     a3 = data.frame(aic = fits.shn$fit[[1]]$aic),
     b1 = fits.shn$info[[2]]$coefficients, 
     b2 = data.frame(r2 = with(fits.shn$info[[2]], 1 - deviance/null.deviance)), 
     b3 = data.frame(aic = fits.shn$fit[[2]]$aic)
) %>% 
    lapply(as.data.frame) %>% 
    lapply(rownames_to_column) %>% 
    writexl::write_xlsx("shn.result.xlsx")

# Ordination --------------------------------------------------------------
dis <-wide %>% 
    mutate(zone = fct_collapse(
        zone, `импактная` = c("импактная", "суперимпактная"))) %>%
    select(-id, -km, -tur) %>% 
    group_by(year, zone, site, plot) %>% 
    summarise_all(mean) %>% 
    ungroup %>% 
    # select_if(~ !is.numeric(.) || sum(.) != 0)
    unite(ID, year, zone, site, plot, sep = "_") %>% 
    filter(rowSums(.[,2:ncol(.)]) > 0 
           # str_detect(ID, "супер", negate = TRUE)
    ) %>%
    # filter(str_detect(ID, "супер", negate = TRUE))
    column_to_rownames("ID") %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
pc <- ape::pcoa(dis)
eig <- pc$values$Eigenvalues
eig <- round(eig/sum(eig)*100, 1)
pc <- pc$vectors %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    select(1:3) %>% 
    separate(ID, into = c("year", "zone", "site", "plot"), 
             sep = "_") %>% 
    as_tibble() %>% 
    mutate(#zone = str_replace_all(zone, "ая", "ый"),
           zone = factor(zone, levels = c("фоновая", "буферная", "импактная")),
           zone = fct_relabel(zone, ~paste0(.x, " территория")))

ggplot(pc, aes(x = Axis.1, y = Axis.2, fill = year, color = year, shape = year)) +
    geom_point(color = "black") +  
    stat_ellipse(mapping = aes(color = year), linewidth = 0.5, linetype = "dashed") +
    facet_grid(cols = vars(zone)) +
    scale_shape_manual(values = c(21, 22))+ 
    scale_fill_manual(values = c("red", "blue")) +
    scale_color_manual(values = c("red", "blue")) +
    labs(subtitle = "По динамической плотности", 
         x = paste0("Ось 1 (", eig[1], " %)"), 
         y = paste0("Ось 2 (", eig[2], " %)"), 
         fill = "Год", color = "Год", shape = "Год") +
    ggplot2::guides(size = "none", color = "none") + 
    theme(panel.grid = element_blank())
ggsave("Plot5.svg", width = 22, height = 9, units = "cm")
# export ------------------------------------------------------------------
list(abundance = fits.abu, 
     number.of.species = fits.nsp, 
     shannon.index = fits.shn) %>% 
    lapply(function(a){select(a, -fit, -info)}) %>% 
    writexl::write_xlsx(paste0("diversity.models_", Sys.Date(), ".xlsx"))

# Rarefication ------------------------------------------------------------
library(parallel)
cl <- makeCluster(detectCores()-1)

rar <- wide %>% 
    unite(id, site, year, zone,  sep = "_") %>% 
    select(-tur, -km, -plot) %>% 
    group_by(id) %>% 
    summarise_all(sum) %>% 
    column_to_rownames("id") %>% 
    t %>% 
    as.data.frame() %>% 
    as.list()
rar <- rar %>% 
    .[which(sapply(rar, function(a){length(a[a>0])})>1)] %>% 
    lapply(function(a){a[a>0]}) %>% 
    parLapply(cl = cl, ., function(a){
        a |> 
            iNEXT::iNEXT(size = seq(0, 200, by = 5), nboot = 0) |>
            purrr::pluck("iNextEst", "size_based") |> 
            dplyr::select(m, Method, qD) |> 
            dplyr::filter(m %in% seq(0, 200, by = 5) | Method == "Observed") 
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    as_tibble() %>% 
    separate(id, into = c("site", "year", "zone"), sep = "_") %>% 
    mutate(zone = factor(zone, ordered = TRUE,
        levels = c("фоновая", "буферная", "импактная", "суперимпактная")))

df <- rar %>% 
    filter(Method == "Observed") %>% 
    group_by(year, zone) %>% 
    summarise(m = mean(m), qD = mean(qD), .groups = "drop") %>% 
    mutate(m = case_when(m > 200 ~ 200, TRUE ~ m))

rar %>% 
    filter(m<=200, Method != "Observed" | zone == "суперимпактная") %>% 
    filter(zone != "суперимпактная" | Method != "Extrapolation") %>% 
    group_by(year, zone, m) %>% 
    summarise(ssd = sd(qD), qD = mean(qD), .groups = "drop") %>% 
    ggplot(aes(x = m, y = qD, 
               color = zone, fill = zone)) + 
    facet_wrap(~year) +
    geom_ribbon(aes(ymin = qD-ssd, ymax = qD+ssd), alpha = 0.2, color = NA) +
    geom_line() + 
    geom_point(mapping = aes(x = m, y = qD,  color = zone, fill = zone),
               data = df, shape = 22) + 
    scale_color_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    scale_fill_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    labs(x = "Individuals", y = "Species") + 
    theme(panel.grid = element_blank())
ggsave("Raref.png", width = 10, height = 5.5, dpi = 600)

# Dominant species --------------------------------------------------------
wide %>% 
    select(id:plot, Pterostichus_oblongopunctatus) %>% 
    # mutate(Pterostichus_oblongopunctatus = log(Pterostichus_oblongopunctatus+1)) %>%
    mutate(Pterostichus_oblongopunctatus = MASS::boxcox(wide$Pterostichus_oblongopunctatus)) %>%
pull(Pterostichus_oblongopunctatus) %>% shapiro.test() 

# Models analysis ---------------------------------------------------------


fits %>% 
    filter(str_detect(ff, "shan", negate = TRUE), 
           # shapiro >=0.05
    ) %>% 
    select(-nm, -fit, -info, -var.res, -var.ftd )







