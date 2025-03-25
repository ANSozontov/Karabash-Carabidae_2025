export <- FALSE
# Data load ---------------------------------------------------------------
library(tidyverse)
theme_set(
    theme_bw() +
        theme(
            text = element_text(family = "serif", size = 15),
            legend.position = "bottom"
        )
)

long <- readxl::read_excel("Data/Carabidae_25.01.2023.xlsx", 
                           sheet = "main_data") %>% 
    filter(taxa != "no_insects") %>% 
    select(-duration, -traps) %>% 
    arrange(desc(site), taxa) %>% # -taxa
    mutate(tur = factor(tur), 
        site = fct_inorder(site),
        zone = factor(zone, 
            levels = c("fon", "bufer", "impact", "superimpact")),
        zone = fct_recode(zone, 
            "background" = "fon",
            "industrial barren" = "superimpact"),
        km = str_extract(site, "[:digit:]{1,}"),
        km = as.numeric(km), 
        year = as.factor(year), 
        .after = "site")

# 1 = turs 1 and 2 are united
wide <- long %>% 
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0)

nsp_100 <- wide %>% 
    select(-km) %>% 
    unite("id", year, zone, site, plot, sep = "_") %>% 
    column_to_rownames("id") %>% 
    t %>% 
    as.data.frame() %>% 
    sapply(function(x){
        x <- x[x>0]
        x <- as.numeric(as.character(x))
        if(length(x) == 0){0} else 
            if(length(x) == 1){1} else {
                iNEXT::iNEXT(x, size = 100, se = FALSE, q = 0) %>%
                    pluck("iNextEst", "size_based") %>%
                    filter(m == 100) %>%
                    pull(qD)
            }
    })

# 1 = turs 1 and 2 are united
div <- tibble(wide[,1:5], 
               abu = apply(wide[,6:ncol(wide)], 1, sum),
               nsp = apply(wide[,6:ncol(wide)], 1, function(a){length(a[a>0])}),
               nsp100 = nsp_100,
               shan= vegan::diversity(wide[,6:ncol(wide)], 
                                      MARGIN = 1, index = "shannon", base = exp(1))
) %>% 
    mutate(km2 = km^2, kmLog = log(km), .after = km) %>% 
    mutate(abuLog = log(abu+1), .after = abu)

rm(nsp_100)
res <- list()
plots <- list(NA)
tables <- list(NA)

# Rarefication ------------------------------------------------------------
library(parallel)
cl <- makeCluster(detectCores()-1)
rar <- long %>% 
    group_by(year, zone, taxa) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    unite("id", year, zone, sep = "_") %>% 
    split(.$id) %>% 
    map(~sort(.x$abu[.x$abu>0])) %>% 
    parLapply(cl = cl, ., function(a){
        a |> 
            as.character() |>
            as.numeric() |>
            iNEXT::iNEXT(
                size = seq(0, 100, by = 2), 
                ### 999
                nboot = 9, 
                ### 999
                se = TRUE, 
                conf = 0.999) |>
            purrr::pluck("iNextEst", "size_based") |>
            dplyr::select(m, Method, qD, qD.LCL,   qD.UCL) 
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    as_tibble() %>% 
    separate(id, into = c("year", "zone"), sep = "_") %>% 
    mutate(zone = factor(zone, ordered = TRUE, levels = levels(long$zone)))

plots$raref <- rar %>% 
    filter(m %in% seq(0, 100, by = 2)) %>% 
    mutate(qD.LCL  = case_when(is.na(qD.LCL) ~ qD, TRUE ~ qD.LCL),
           qD.UCL  = case_when(is.na(qD.UCL) ~ qD, TRUE ~ qD.UCL)) %>% 
    ggplot(aes(x = m, y = qD, color = zone, fill = zone)) + 
    facet_wrap(~year) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.2, color = NA) +
    geom_line() + 
    scale_color_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    scale_fill_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    labs(x = "Individuals", y = "Species") + 
    theme(panel.grid = element_blank())

if(!export){plots$raref}

# Models template ---------------------------------------------------------
models_fit <- function(x){
    x %>% 
        split(1:nrow(.)) %>%
        lapply(function(y){
            df <-  filter(div, year == y$year)
            if(str_detect(y$formula, "Segmented")){
                ###
                segmented::segmented(lm(
                    str_replace_all(y$formula, "Segmented", ""), 
                    data = df), seg.Z = ~km)
                ###
            } else {
                lm(formula = y$formula, data = df)
            }
        })
}

models_pred <- function(fits){
    lapply(fits, function(x){
        c(if("segmented" %in% class(x)){x$psi[2]}else{numeric()}, seq(1, 32, by = 0.5)) %>% 
            sort() %>% 
            unique() %>%
            tibble(
                km = .,
                km2 = .^2, 
                kmLog = log(.)
            ) %>% 
            mutate(predicted = predict(x, .))
    }
)}
    
        mutate(
            linear    = predict(fits$fit[[1]], .),
            segmented = predict(fits$fit[[2]], .),
            nonlinear = predict(fits$fit[[3]], .), 
            logarithm = predict(fits$fit[[4]], .)
        ) %>% 
        dplyr::select(-km2, -kmLog)
    if(str_detect(toupper(formulas[1]), "LOG")) {
        pred <- pred %>% 
            mutate_at(3:ncol(.), exp)
    }
    pred <- pred %>% 
        `colnames<-`(c("year", "km", formulas)) %>% 
        pivot_longer(names_to = "model", values_to = "abu", -c("year", "km"))
    pred$model <-  map_chr(str_split(pred$model, "year \\+ "), ~.x[2])
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

model_viz <- function(df, yy){
    div1 %>% 
        rename(yy = which(colnames(div1)==yy)) %>% 
        ggplot(aes(x = km, y = yy)) + 
        geom_line(aes(km, abu, color = model), data = map_dfr(df$d, rbind), 
                  linewidth = 1) + 
        geom_point(shape = 21, size = 3) + 
        facet_wrap(~year, scales = "fixed", ncol = 2) +
        theme(panel.grid = element_blank()) +
        guides(fill="none")
}

# Abundance ---------------------------------------------------------------
# res$abundance <- 
expand_grid(
        formula = paste0("abu ~ ", c("km", "kmSegmented", "km + km2", "kmLog")), 
        year = c(2009, 2014)) %>% 
    mutate(
        fits = models_fit(.), 
        pred = models_pred(fits), 
        aic = map_dbl(fits, ~AIC(.x)), 
        r2 = map_dbl(fits, ~summary(.x)$adj.r.squared)
           ) # %>% arrange(year, aic)

res$abundance
model_viz(res$abundance, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели исходных показателей обилия")

# ggsave("1a. Abundance.png", height = 8, width = 11, dpi = 600)

# Abundance LLOG -----------------------------------------------------------
res$abundance_log <- c("km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("abuLog ~ year + ", .) %>% 
    models_fit()
res$abundance_log
model_viz(res$abundance_log, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели логарифмированных показателей обилия")
# ggsave("1b. Abundance_log.png", height = 8, width = 11, dpi = 600)

# N_species ---------------------------------------------------------------
res$nsp <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4", 
    paste0("nsp ~ year + ", .) %>% 
    models_fit()
res$nsp
model_viz(res$nsp, "nsp") + 
    labs(x = NULL, y = "Количество видов",
         subtitle = "Видовое богатство")

# ggsave("2. n_species.png", height = 8, width = 11, dpi = 600)

# N_species rarefication --------------------------------------------------
res$nsp100 <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4", 
    paste0("nsp100 ~ year + ", .) %>% 
    models_fit()
res$nsp100
model_viz(res$nsp100, "nsp100") + 
    labs(x = NULL, y = "Количество видов",
         subtitle = "Видовое богатство (разрежение: 100 экз.)")

# Shannon -----------------------------------------------------------------
res$shan <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4",
    paste0("shan ~ year + ", .) %>% 
    models_fit()
res$shan
model_viz(res$shan, "shan") + 
    labs(x = NULL, y = "Индекс Шеннона",
         subtitle = "Видовое разнообразие")
# ggsave("3. Shannon.png", height = 8, width = 11, dpi = 600)

# final export ------------------------------------------------------------
# rarefaction fig
ggsave(paste0("export/Fig.x_raref_", Sys.Date(), ".pdf"), plots$raref, 
       width = 9, height = 5.5, dpi = 600)


# Fig. 2
div1 %>% 
    select(year, km, 
           A_abuLog = abu, C_nsp = nsp, 
           D_nsp100 = nsp100, B_shan = shan) %>% 
    pivot_longer(names_to = "type", values_to = "abu", -1:-2) %>% 
    ggplot(aes(km, abu, color = year)) + 
    geom_line(
        # linetype = "dashed",
        data = mutate(
            rbind(
                res$abundance_log$d[[2]], res$shan$d[[2]],
                res$nsp$d[[2]], res$nsp100$d[[2]]),
            type = rep(c("A_abuLog", "B_shan", "C_nsp", "D_nsp100"), each = 128))
    ) +
    geom_point(shape = 21, size = 2) +
    facet_wrap(
        ~type,
        scales = "free") + 
    labs(x = "Distance, km", y = NULL, color = "Year")
ggsave(paste0("export/Fig.2_segm_", Sys.Date(), ".pdf"), 
       width = 6.5, height = 5.5, dpi = 600)

# Tables 
res %>% 
    map(~select(.x, -d, -fit)) %>% 
    writexl::write_xlsx(
        paste0("export/models_all_", Sys.Date(), ".xlsx")
    )

all_fits <- res %>% 
    map(~dplyr::select(.x[2,], fit)[[1]]) %>% 
    map(~.x[[1]])

all_fits %>% 
    lapply(summary) %>% 
    lapply(capture.output) %>% 
    map(~tibble(results = .)) %>% 
    map2(
        .,
        all_fits %>% 
            lapply(function(segmented_model){
                
                # psi
                est <- segmented_model$psi["psi1.km", "Est."]
                se <- segmented_model$psi["psi1.km", "St.Err"]
                t.val <- est / se
                psi_p.val <- 2 * pt(-abs(t.val), df = df.residual(segmented_model))
                if(psi_p.val<0.0001){
                    psi_p.val <- "<0.0001"
                } else {
                    psi_p.val <- round(psi_p.val)
                }
                
                #km
                coefs <- summary(segmented_model)$coefficients
                covariance <- vcov(segmented_model)["km", "U1.km"]
                
                coef2 <- coefs["km","Estimate"] + coefs["U1.km","Estimate"]
                se_coef2 <- sqrt(coefs["km","Std. Error"]^2 + coefs["U1.km","Std. Error"]^2 + 2 * covariance)
                t.val <- coef2 / se_coef2
                U1.km_p.val <- 2 * pt(-abs(t.val), df = df.residual(segmented_model) )
                if(U1.km_p.val<0.0001){
                    U1.km_p.val <- "<0.0001"
                } else {
                    U1.km_p.val <- round(U1.km_p.val, 4)
                }
                
                #return
                return(tibble(results = c(
                    rep("", 2),
                    "Custom results:",
                    paste0("psi1.km p-value = ", psi_p.val),
                    paste0("U1.km coeff = ", round(coef2, 3)),
                    paste0("U1.km coeff SE = ", round(se_coef2, 3)),
                    paste0("U1.km coeff_p.val = ", U1.km_p.val)
                    )))
            }),
        ~rbind(.x, .y)) %>% 
    writexl::write_xlsx(paste0("export/models_selected_", Sys.Date(), ".xlsx"))


# Supplement 3 ------------------------------------------------------------
p <- gridExtra::grid.arrange(
    model_viz(res$abundance_log, "abu") + 
        labs(x = NULL, y = NULL, #"Обилие (особей на 100 лов.-сут.)", 
             subtitle = "1") + 
        theme(legend.position = "none"),
    
    model_viz(res$nsp, "nsp") + 
        labs(x = NULL, y = NULL, #"Количество видов",
             subtitle = "2") + #Видовое богатство"
        theme(legend.position = "none"),
    
    model_viz(res$nsp100, "nsp100") + 
        labs(x = NULL, y = NULL, #"Количество видов",
             subtitle = "3") + #"Видовое богатство (разрежение: 100 экз.)")+ 
        theme(legend.position = "none"),
    
    model_viz(res$shan, "shan") + 
        scale_color_discrete(labels = LETTERS[1:4]) + 
        labs(x = NULL, y = NULL, #"Индекс Шеннона",
             color = "Model", 
             subtitle = "4"), #Видовое разнообразие"),
    ncol = 1, widths = c(1)
)

ggsave(paste0("export/Suppl.3_", Sys.Date(), ".pdf"), plot = p, width = 8, height = 12)

# Multidimensional count --------------------------------------------------
dis <- wide %>% 
    select(-site) %>% 
    unite("ID", zone, year, km, plot, sep = "_") %>% 
    column_to_rownames("ID") %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
pc <- ape::pcoa(dis)
eig <- pc$values$Eigenvalues
eig <- round(eig/sum(eig)*100, 1)
pc <- pc$vectors %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    as_tibble() %>% 
    separate(ID, into = c("zone", "year", "site", "plot"), 
             sep = "_") %>% 
    mutate(zone = factor(zone, levels = levels(long$zone)))

L <- expand_grid(levels(long$zone), c(2009, 2014)) %>% 
    apply(1, function(a){paste0(a, collapse = "_")})

distances <- as.matrix(dis)
distances[upper.tri(distances)] <- NA
diag(distances) <- NA
distances <- distances %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    as_tibble() %>%
    pivot_longer(names_to = "id2", values_to = "dis", -id1) %>%
    filter(!is.na(dis)) %>%
    separate(id1, c("zone1", "year1"), sep = "_", extra = "drop") %>%
    separate(id2, c("zone2", "year2"), sep = "_", extra = "drop") %>%
    mutate(
        id1 = factor(paste0(zone1, "_", year1), levels = L, ordered = TRUE),
        id2 = factor(paste0(zone2, "_", year2), levels = L, ordered = TRUE),
        .keep = "unused") %>% 
    split(1:nrow(.)) %>% 
    lapply(function(x){
        if(x$id1 > x$id2) {
            tibble(dis = x$dis, id1 = x$id2, id2 = x$id1)
        } else {
            tibble(dis = x$dis, id1 = x$id1, id2 = x$id2)
        }
    }) %>% 
    map_dfr(rbind) %>% 
    pivot_wider(
        names_from = id2, values_from = dis, 
        values_fill = NA, values_fn = mean)
distances %>% 
    mutate_if(is.numeric, function(a){round(a, 2)}) %>% 
    writexl::write_xlsx(paste0("export/distances_", Sys.Date(), ".xlsx"))

# Multidimensional viz ----------------------------------------------------
ggplot(pc, aes(x = Axis.1, y = Axis.2, linetype = year,
        fill = zone, color = zone, shape = year)) +
    geom_point(color = "black", size = 2.5) +  
    stat_ellipse() + 
    scale_shape_manual(values = c(21, 22))+
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_fill_manual(values = c("darkgreen", "greenyellow", "orange", "red")) +
    scale_color_manual(values = c("darkgreen", "greenyellow", "orange", "red")) +
    labs(subtitle = "Ординация по динамической плотности\nТуры объединены", 
         x = paste0("Ось 1 (", eig[1], " %)"), 
         y = paste0("Ось 2 (", eig[2], " %)"), 
         # fill = "Год", color = "Год", shape = "Год"
    ) +
    theme(panel.grid = element_blank())

ggsave(paste0("export/Fig.3_ord_", Sys.Date(), ".svg"), width = 18, height = 13, units = "cm")
ggsave(paste0("export/Fig.3_ord_", Sys.Date(), ".png"), width = 18, height = 13, units = "cm")

