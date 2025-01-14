# data load ---------------------------------------------------------------
res <- list()
library(tidyverse)
theme_set(
    theme_bw() + 
    theme(
        legend.position = "bottom", 
        text = element_text(family = "serif", size = 20),
        strip.background = element_rect(fill = "aliceblue")
    ))
L <- rio::import("https://docs.google.com/spreadsheets/d/1KmTMO-wg95E7u8faRWDfaze3_PenEKCtMMtVZvLx6hw/edit?gid=1059328530") %>% 
    as_tibble() %>% 
    select(year:mocerat) %>% 
    pivot_longer(names_to = "sexage", values_to = "abu", -year:-species) %>% 
    filter(abu>0) %>% 
    mutate(
        tur = tur-1,
        km = as.numeric(str_extract(site, "[:digit:]+")),
        species = str_squish(species),
        taxa = case_when(substr(species, nchar(species)-2, nchar(species)) == " sp" 
            ~ paste0(species, " (", family, ")"), 
            TRUE ~ species),
        abu = abu/traps/5*100, 
        .after = site) %>% 
    transmute(year = as.factor(year), tur, km, plot, defects, taxa, sexage, abu)
    
div <- L %>% 
    unite("id", year, tur, plot, km, sep = "_") %>% 
    split(.$id)
   

div %>% map_dbl(~.x %>% pull(abu) %>% sum) %>% tibble(x = .) %>% ggplot(aes(x)) + geom_density()

div <- div %>% 
    # `[`(100:118) %>%
    lapply(function(x){
        abu <-  x %>% 
            filter(defects == 0) %>% 
            pull(abu) %>% 
            sum()
        if(nrow(x) < 2) {
            return(tibble(
                abu = abu, 
                nsp  = 1, 
                nsp50 = NA, 
                shan = 0, 
                iBP = 1
            ))
        }
        
        taxa.enabled <- x$taxa %>% 
            unique() %>% 
            sort %>%
            str_split_fixed(" ", 2) %>%
            `colnames<-`(c("gen", "sp")) %>% 
            as_tibble() %>% 
            split(substr(.$sp, 1, 3) == "sp ") %>% 
            `names<-`(c("sp", "gen")[1:length(.)])
        
        if(length(taxa.enabled) == 1) {
            taxa.enabled <- taxa.enabled[[1]] %>% 
                unite("sp", gen, sp, sep = " ") %>% 
                pull(sp)
        } else { 
            taxa.enabled <- taxa.enabled$gen %>% 
                filter(!(gen %in% taxa.enabled$sp$gen)) %>% 
                rbind(taxa.enabled$sp, .) %>% 
                unite("sp", gen, sp, sep = " ") %>% 
                pull(sp)
        }
        
        x <- x %>% 
            filter(taxa %in% taxa.enabled, defects == 0, abu > 0) %>% 
            group_by(taxa) %>% 
            summarise(abu = sum(abu))
        
        tibble(
            abu = abu, 
            nsp  = length(taxa.enabled), # sort() %>% paste0(collapse = ",")
            nsp50 = x %>% 
                pull(abu) %>% 
                ceiling %>% 
                iNEXT::iNEXT(q = 0, size = 100, se = FALSE, nboot = 0)  %>% 
                pluck("iNextEst", "size_based") %>% 
                filter(m == 100) %>% 
                pull(qD), 
            shan = x %>% 
                pull(abu) %>% 
                vegan::diversity(index = "shannon"), 
            iBP = x %>% 
                mutate(abu = abu/sum(abu)) %>% 
                filter(abu == max(abu)) %>% 
                pull(abu)
            )
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    distinct() %>% 
    separate(id, into = c("year", "tur", "plot", "km")) %>% 
    mutate_at(1:4, as.numeric) %>% 
    mutate(km2 = km^2, kmLog = log(km), year = factor(year), tur = factor(tur), abuLog = log(abu+1), .after = km)

# Models template ---------------------------------------------------------
models_fit <- function(formulas){
    fits <- tibble(ff = formulas) %>% 
        split(1:nrow(.)) %>% 
        lapply(function(a){ 
            if(str_detect(a$ff, "Segmented")){
                #######
                segmented::segmented(lm(
                    str_replace_all(a$ff, "Segmented", ""), 
                    data = div), seg.Z = ~km)
                #######
            } else {
                lm(formula = a$ff, data = div)
            }
        }) %>% 
        tibble(ff = formulas, fit = .)
    pred <- seq(1, 32, by = 0.5) %>% 
        c(fits$fit[[2]]$psi[2]) %>% 
        sort() %>% 
        unique() %>%
        expand_grid(
            year = factor(c(2009, 2014)),
            tur = factor(c(1, 2)),
            km = .) %>%
        mutate(
            km2 = km^2, 
            kmLog = log(km)
        ) %>% 
        mutate(
            linear    = predict(fits$fit[[1]], .),
            segmented = predict(fits$fit[[2]], .),
            nonlinear = predict(fits$fit[[3]], .), 
            logarithm = predict(fits$fit[[4]], .)
        ) %>% 
        select(-km2, -kmLog)
    if(str_detect(toupper(formulas[1]), "LOG")) {
        pred <- pred %>% 
            mutate_at(4:ncol(.), exp)
    }
    pred <- pred %>% 
        `colnames<-`(c("year", "tur", "km", formulas)) %>% 
        pivot_longer(names_to = "model", values_to = "abu", -c("year", "tur", "km"))
    pred$model <-  map_chr(str_split(pred$model, "year \\+ tur \\+ "), ~.x[2])
    pred <- pred %>% 
        mutate(model = factor(model, 
                              levels = c("km", "kmSegmented", "km + km2", "kmLog")) # "zone3", "zone4", 
        ) %>% 
        split(.$model)
    
    fits %>% 
        pull(fit) %>% 
        lapply(function(a){
            tibble(
                aic = AIC(a), 
                r2 = summary(a)$adj.r.squared, 
                sh.test = round(shapiro.test(a$residuals)$p.value, 4)
            ) 
        }) %>% 
        map_df(rbind) %>% 
        tibble(fits, .) %>% 
        mutate(d = pred)
}

# VIZ template ------------------------------------------------------------
model_viz <- function(df, yy){
    div %>% 
        rename(yy = which(colnames(div)==yy)) %>% 
        ggplot(aes(x = km, y = yy)) + 
        geom_line(aes(km, abu, color = model), data = map_dfr(df$d, rbind),
                  linewidth = 1, linetype = "solid", alpha = 0.7) +
        geom_point(shape = 21, size = 3) + 
        facet_grid(cols = vars(year), rows = vars(tur)) + 
        theme(panel.grid = element_blank()) +
        guides(fill="none")
}

# Models ------------------------------------------------------------------
# Abundance
res$abundance <- c("km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("abu ~ year + tur + ", .) %>% 
    models_fit()
res$abundance
model_viz(res$abundance, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели исходных показателей обилия\n", 
         color = "Тип модели")
ggsave("export/1a. Abundance.png", height = 8, width = 11, dpi = 600)

# Abundance LOG 
res$abundance_log <- c("km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("abuLog ~ year + tur + ", .) %>% 
    models_fit()
res$abundance_log
model_viz(res$abundance_log, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели логарифмированных показателей обилия\n", 
         color = "Тип модели")
         
ggsave("export/1b. Abundance_log.png", height = 8, width = 11, dpi = 600)

# N_species 
res$nsp <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4", 
    paste0("nsp ~ year + tur + ", .) %>% 
    models_fit()
res$nsp
model_viz(res$nsp, "nsp") + 
    labs(x = NULL, y = "Количество видов",
         subtitle = "Видовое богатство\n", 
         color = "Тип модели")
ggsave("export/2. n_species.png", height = 8, width = 11, dpi = 600)

# N_species d50  interpolated - extrapolated 
res$nsp50 <- c("km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("nsp50 ~ year + tur + ", .) %>% 
    models_fit()
res$nsp50
model_viz(res$nsp50, "nsp") + 
    labs(x = NULL, y = "Количество видов на 50 особей",
         subtitle = "Видовое богатство (на 50 особей)\n", 
         color = "Тип модели")
ggsave("export/3. n_species rar-ext to 50.png", height = 8, width = 11, dpi = 600)

# Shannon 
res$shan <- c("km", "kmSegmented", "km + km2", "kmLog") %>%
    paste0("shan ~ year + tur + ", .) %>% 
    models_fit()
res$shan
model_viz(res$shan, "shan") + 
    labs(x = NULL, y = "Индекс Шеннона",
         subtitle = "Видовое разнообразие\n", 
         color = "Тип модели")
ggsave("export/4. Shannon.png", height = 8, width = 11, dpi = 600)

# Berger-Parker 
# res$iBP <- c("km", "kmSegmented", "km + km2", "kmLog") %>% 
#     paste0("iBP ~ year + tur + ", .) %>% 
#     models_fit()
# res$iBP
# model_viz(res$iBP, "iBP") + 
#     labs(x = NULL, y = "Индекс Бергера-паркера",
#          subtitle = "Доминирование")
# ggsave("export/5. Berger-Parker.png", height = 8, width = 11, dpi = 600)

# Multidimensional --------------------------------------------------------
pcoa_res <- list()
pcoa_custom <- function(w, binary.form = FALSE){
    if("zone" %in% colnames(w)) {
        w <- select(w, -zone)
    }
    dis <- vegan::vegdist(w, method = "bray", binary = binary.form)
    pcoa <- ape::pcoa(dis)
    eig <- pcoa$values$Eigenvalues
    if(min(eig) < 0){
        eig <- eig + abs(min(eig))
    }
    eig <- round(eig/sum(eig)*100, 1)
    pcoa <- pcoa$vectors %>%
        as_tibble() %>%
        transmute(id = row.names(pcoa$vectors), Axis.1, Axis.2) %>%
        separate(id, into = c("year", "tur", "plot", "km"), sep = "_") %>%
        mutate_all(as.numeric) %>% 
        mutate(zone = case_when(km > 25 ~ "А. Фоновая", km >= 9 ~ "Б. Буферная", TRUE ~ "В. Импактная"),
               .after = km)
    lst(dis, pcoa, eig)
}

pcoa_viz <- function(x, st = "", xx = 1, yy = 1){ 
    x$pcoa %>% 
    mutate(
        year = as.factor(year),
        tur = as.factor(tur),
        Axis.1 = Axis.1 * xx, 
        Axis.2 = Axis.2 * yy,
        ) %>% 
    ggplot(aes(x = Axis.1, y = Axis.2, color = zone, shape = tur)) + 
        geom_point(size = 3, alpha = 0.6)+
        stat_ellipse(level = 0.9) +
        scale_color_manual(values = c("green4", "gold2", "red")) + 
        scale_shape_manual(values = c(21, 25)) + 
        labs(subtitle = st, 
             x = paste0("Ось 1, ", x$eig[1], "%"),
             y = paste0("Ось 2, ", x$eig[2], "%"), 
             color = "Зона", 
             shape = "Тур"
             ) + 
        facet_grid(cols = vars(year)) #, rows = vars(tur)) 
    }

wide <- L %>%
    filter(defects == 0) %>%
    select(-defects, -sexage) %>%
    # mutate(zone = case_when(km > 25 ~ "Фоновая", km >= 9 ~ "Буферная", TRUE ~ "Импактная")) %>% 
    unite("id", year, tur, plot, km, sep = "_") %>%
    arrange(taxa, id) %>% 
    pivot_wider(names_from = taxa, values_from = abu, values_fill = 0, values_fn = sum) %>%
    column_to_rownames("id")

# raw data
pcoa_res$raw <- wide %>% pcoa_custom()
# pcoa_res$raw$pcoa %>% pull(Axis.2) %>% range
pcoa_viz(pcoa_res$raw, "Население. Обилия\n")
    # scale_x_continuous(limits = c(-0.44, 0.37)) +
    # scale_y_continuous(limits = c(-0.3, 0.54))
ggsave("export/A. Обилия.png", height = 8, width = 11, dpi = 600)
pcoa_res$raw$permanova <- vegan::adonis2(
    pcoa_res$raw$dis ~ tur + year + zone, 
    data = pcoa_res$raw$pcoa, 
    permutations = 9999,
    by = "terms")

# freq data
pcoa_res$freq <- wide %>% 
    t %>% 
    as.data.frame() %>% 
    mutate_all(function(a){a/sum(a)}) %>% 
    t %>% 
    as.data.frame() %>% 
    pcoa_custom()
# pcoa_res$freq$pcoa %>% pull(Axis.2) %>% range
pcoa_viz(pcoa_res$freq, "Население. Доли\n") 
    # scale_x_continuous(limits = c(-0.44, 0.4)) +
    # scale_y_continuous(limits = c(-0.3, 0.54))
ggsave("export/B. Доли.png", height = 8, width = 11, dpi = 600)
pcoa_res$freq$permanova <- vegan::adonis2(
    pcoa_res$freq$dis ~ tur + year + zone, 
    data = pcoa_res$freq$pcoa, 
    permutations = 9999,
    by = "terms")

# binary
pcoa_res$bin <- wide %>% pcoa_custom(binary.form = TRUE)
# pcoa_res$bin$pcoa %>% pull(Axis.2) %>% range
pcoa_viz(pcoa_res$bin, "Состав (бинарный)\n", 1, 1)
    # scale_x_continuous(limits = c(-0.31, 0.5)) +
    # scale_y_continuous(limits = c(-0.33, 0.54))
ggsave("export/C. Состав.png", height = 8, width = 11, dpi = 600)
pcoa_res$bin$permanova <- vegan::adonis2(
    pcoa_res$bin$dis ~ tur + year + zone, 
    data = pcoa_res$bin$pcoa, 
    permutations = 9999,
    by = "terms")

# export
pcoa_res %>% 
    map(~pluck(.x, "permanova")) %>% 
    lapply(capture.output) %>% 
    map(~tibble(results = .)) %>% 
    writexl::write_xlsx(paste0("export/permanova_", Sys.Date(), ".xlsx"))

# final export ------------------------------------------------------------
spider_dictionary <- L %>% 
    select(taxa) %>% 
    distinct() %>% 
    filter(str_detect(taxa, " sp ", negate = TRUE), 
           str_detect(taxa, "[:digit:]", negate = TRUE), 
           !(taxa %in% c("Lacinius ephippiatus",
                       "Lophopilio palpinalis", 
                       "Nemastoma lugubre", 
                       "Oligolophus tridens"))
           ) %>% 
    arrange(taxa) %>% 
    split(1:nrow(.)) %>% 
    map(~rgbif::name_backbone(.x$taxa[1], rank = "SPECIES", kingdom = "Animalia", order = "Araneae")) %>% 
    map_dfr(rbind)
    
spider_dictionary <- spider_dictionary %>% 
    select(taxa = canonicalName, taxa_a = scientificName )

res$abundance_d100 <- L %>% 
    left_join(spider_dictionary, by = "taxa") %>% 
    mutate(
        taxa = case_when(is.na(taxa_a) ~ taxa, TRUE ~ taxa_a),
        zone = case_when(km > 25 ~ "А. Фоновая", km >= 9 ~ "Б. Буферная", TRUE ~ "В. Импактная"),
        .after = km) %>% 
    group_by(year, zone, taxa) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    mutate(abu = round(abu)) %>% 
    pivot_wider(names_from = year, values_from = abu) %>% 
    mutate_at(3:4, function(a){a <- as.character(a); a[is.na(a)] <- "-"; a}) %>% 
    mutate(abu = paste0(`2009`, "|", `2014`), .keep = "unused") %>% 
    pivot_wider(names_from = zone, values_from = abu, values_fill = "-") 

res$abundance_raw <- 
    rio::import("https://docs.google.com/spreadsheets/d/1KmTMO-wg95E7u8faRWDfaze3_PenEKCtMMtVZvLx6hw/edit?gid=1059328530") %>% 
    as_tibble() %>% 
    mutate_at(c("m", "f", "sm", "sf", "j", "mocerat"), function(a){a[is.na(a)] <- 0;a}) %>% 
    transmute(
        year, 
        taxa = species, 
        abu = m + f + sm + sf + j + mocerat, 
        zone = as.numeric(str_extract(site, "[:digit:]+")), 
        zone = case_when(zone > 25 ~ "А. Фоновая", zone >= 9 ~ "Б. Буферная", TRUE ~ "В. Импактная")
    ) %>% 
    left_join(spider_dictionary, by = "taxa") %>% 
    mutate(taxa = case_when(is.na(taxa_a) ~ taxa, TRUE ~ taxa_a)) %>% 
    group_by(year, zone, taxa) %>% 
    summarise(abu = sum(abu, na.rm = TRUE), .groups = "drop") %>% 
    unite("i", year, zone, sep = ". ") %>% 
    mutate(abu = as.character(abu)) %>% 
    pivot_wider(names_from = i, values_from = abu, values_fill = "-")
    # mutate_at(3:4, function(a){a <- as.character(a); a[is.na(a)] <- "-"; a}) 
    # mutate(abu = paste0(`2009`, "|", `2014`), .keep = "unused") %>% 
    # pivot_wider(names_from = zone, values_from = abu, values_fill = "-") %>% 
    # arrange(taxa)

writexl::write_xlsx(res$abundance_raw, "export/supplement_table.xlsx")
    
    
    
    


res %>% 
    `[`(1:5) %>% 
    map(~select(.x, -d, -fit)) %>% 
    map(~mutate(.x, aic = round(aic, 1))) %>% 
    map(~mutate(.x, r2 = round(r2, 2))) %>% 
    writexl::write_xlsx(paste0("export/models_all_", Sys.Date(), ".xlsx"))

res %>% 
    `[`(1:5) %>% 
    map(~select(.x[2,], fit)[[1]]) %>% 
    map(~.x[[1]]) %>% 
    lapply(summary) %>% 
    lapply(capture.output) %>% 
    map(~tibble(results = .)) %>% 
    writexl::write_xlsx(paste0("export/models_selected_", Sys.Date(), ".xlsx"))


