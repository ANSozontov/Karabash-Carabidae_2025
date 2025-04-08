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
    select(-duration, -traps, -num) %>% 
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
    select(-tur) %>% 
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
    mutate(abuLog = log10(abu+1), .after = abu)

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

# old manual AIC
# manual.aic <- function(a){
#     if(!(a$family$family %in% c("quasipoisson", "poisson"))){return(NA)}
#     loglik <- sum(dpois(a$model[,1], fitted.values(a), log = TRUE))
#     phi <- summary(a)$dispersion
#     -2 * loglik + 2 * summary(a)$df[3] * phi
# }

# derived from "AIC.default"
manual_aic <- function (a, df_observed = NA, k = 2, simple = TRUE){
    lls <- logLik(a)
    # lls is logLik observer
    if(is.na(df_observed)){
        df_observed <- attr(lls, "df")
    }
    if(simple){
        -2 * as.numeric(lls) + k * df_observed
    } else {
        data.frame(
            df = df_observed, 
            AIC = -2 * as.numeric(lls) + k * df_observed)
    }
}

models_coeffs <- function(a, multip = 2){
b <- a %>% 
    filter(str_detect(formula, "Segment")) %>% 
    select(year, fits)
b <- b %>% 
    pluck("fits") %>% 
    `names<-`(b$year) %>% 
    lapply(summary) %>% 
    map(~list(.x$coefficients, .x$psi)) %>% 
    map(~map(.x, ~rownames_to_column(as.data.frame(.x), "pred"))) %>% 
    map_dfr(~rbind(
            select(.x[[1]], pred, `Est.` = 2, `St.Err` = 3), 
            transmute(.x[[2]], pred, `Est.`, `St.Err`)
            ) %>% 
            filter(Est. != 0, pred %in% c("km", "psi1.km")), 
        .id = "year"
        ) %>% 
    split(.$pred) %>% 
    map(~arrange(.x, Est.))

b %>% 
    map_dfr(rbind) %>% 
    mutate_if(is.numeric, ~round(.x, 3)) %>% 
    mutate(CI = paste0("(", Est. - multip*St.Err, "; ", Est. + multip*St.Err, ")"),
           E = paste0(Est., "±", St.Err), 
           .keep = "unused") 
#     
#     map(~.x %>% 
#             mutate_if(is.numeric, ~round(.x, 3)) %>% 
#             unite("i", Est., St.Err, sep = "±") %>% 
#             pivot_wider(names_from = year, values_from = i)
#     )
# r2 <- b %>% 
#     map(~ .x$Est.[2] - multip*.x$St.Err[2] - .x$Est.[1] - multip*.x$St.Err[1]) %>% 
#     map(~tibble(dd = .x))
# 
# list(r1, r2) %>% 
#     transpose() %>% 
#     map_dfr(~cbind(.x[[1]], .x[[2]]))
}

model_viz_supp <- function(df0){
    df1 <- df0 %>% 
        unite("id", formula, year, sep = " ~ ") %>% 
        pull(id) %>% 
        `names<-`(df0$pred, .) %>% 
        map(~select(.x, km, predicted)) %>% 
        map_dfr(rbind, .id = "model") %>% 
        separate(model, into = c("response", "model", "year"), sep = " ~ ")
    
    df2 <- div %>% 
        select_at(c("km", predicted = df1$response[1], "year")) %>% 
        mutate(model = NA)
    
    df1 %>% 
        ggplot(aes(x = km, y = predicted, color = model, linetype = model)) +
        geom_point(shape = 21, size = 3, data = df2, color = "black", alpha = 0.5) +
        geom_line() + 
        facet_wrap(~year, scales = "fixed", ncol = 2)  +
        theme(panel.grid = element_blank()) +
        guides(linetype = "none")
}

model_viz_text <- function(df0){
    df0 <- filter(df0, str_detect(formula, "Segmented")) 
    df1 <- df0 %>% 
        unite("id", formula, year, sep = " ~ ") %>% 
        pull(id) %>% 
        `names<-`(df0$pred, .) %>% 
        map(~select(.x, km, predicted)) %>% 
        map_dfr(rbind, .id = "model") %>% 
        separate(model, into = c("response", "model", "year"), sep = " ~ ")
    
    df2 <- div %>% 
        select_at(c("km", predicted = df1$response[1], "year"))
    
    df3 <- transmute(
        df0, 
        year = as.character(year), 
        xi = map_dbl(df0$fits, ~summary(.x)$psi[2]), 
        ci = map_dbl(df0$fits, ~summary(.x)$psi[3])*2)
    
    df1 %>% 
        ggplot(aes(x = km, y = predicted, color = year, shape = year)) + # , linetype = year
        geom_vline(
            mapping = aes(xintercept = xi, color =  year), 
            data = df3, alpha = 0.8, linetype = "dashed") +
        geom_vline(
            mapping = aes(xintercept = xi+ci, color =  year), 
            data = df3, alpha = 0.8, linetype = "dotted") +
        geom_vline(
            mapping = aes(xintercept = xi-ci, color =  year), 
            data = df3, alpha = 0.8, linetype = "dotted") +
        geom_point(size = 3, data = df2, alpha = 0.5) +
        geom_line() + 
        scale_x_continuous(limits = c(0.1, 32)) +
        theme(panel.grid = element_blank()) +
        guides(linetype = "none")
}

# Models calculation ------------------------------------------------------
models <- expand_grid(
    resp = c("abu", "abuLog", "nsp", "nsp100", "shan"),
    formula = c("km", "kmSegmented", "km + km2", "kmLog"), 
    year = c(2009, 2014)) %>% 
    mutate(
        formula = paste0(resp, " ~ ", formula)) %>% 
    mutate(fits = models_fit(.), 
        pred = models_pred(fits),
        aic = map_dbl(fits, ~AIC(.x)), 
        r2 = map_dbl(fits, ~summary(.x)$adj.r.squared)
    ) %>% 
    mutate(aic2 = map_dbl(fits, ~manual_aic(.x, df_observed = 3)), 
           aic2 = case_when(str_detect(formula, "Segm") ~ aic2, TRUE ~ aic),
           .after = aic) %>% 
    arrange(resp, year, aic2) %>% 
    split(.$resp)

if(!export){models}
labs <- list(
    labs(x = NULL, y = "individuals per 100 traps-days", 
         subtitle = "Abundance"), 
    labs(x = NULL, y = "Log10 of individuals per 100 traps-days", 
         subtitle = "log Abundance"), 
    labs(x = NULL, y = "Number of species",
         subtitle = "Number of species"), 
    labs(x = NULL, y = "Number of species per 100 individuals",
         subtitle = "Number of species rarefied to 100 individuals"),
    labs(x = NULL, y = "Shannon index", subtitle = "Diversity")  
)
plots$suppl <- map2(models, labs, ~model_viz_supp(.x)+.y)

if(!export){plots$suppl}

plots$text <- map2(models, labs, ~model_viz_text(.x)+.y)

if(!export){plots$text}

tables$comps <- map(models, ~models_coeffs(.x, 2))
if(!export){tables$comps}

tables$comps_wide <- tables$comps %>% 
    map(~.x %>% 
            pivot_longer(names_to = "i", values_to = "val", -1:-2) %>% 
            unite("pred", pred, i, sep = "_") %>% 
            pivot_wider(names_from = pred, values_from = val)
    )
if(!export){tables$comps_wide}
            

# Supplement 3 ------------------------------------------------------------

gridExtra::grid.arrange(
    plots$abundance + guides(color = "none"), 
    plots$abundance_log + guides(color = "none"), 
    plots$nsp + guides(color = "none"), 
    plots$nsp100 + guides(color = "none"), 
    plots$shan + guides(color = "none"), 
    ncol = 1
)

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

# Multidimensional  --------------------------------------------------
# Multidimensional count
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

# Multidimensional viz 
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


# Export ------------------------------------------------------------

if(export){
    # rarefaction fig
    ggsave(
        paste0("export/Fig.x_raref_", Sys.Date(), ".pdf"), 
        plots$raref, 
        width = 9, height = 5.5, dpi = 600)
    
    # models tables
    models %>% 
        map(~select(.x, -resp, -fits, -pred)) %>% 
        map(~mutate_if(.x, is.numeric, ~round(.x, 2))) %>% 
        writexl::write_xlsx(
            paste0("export/models_all_", Sys.Date(), ".xlsx")
        )
    
    tables$comps %>% 
        writexl::write_xlsx(
            paste0("export/models_comp_", Sys.Date(), ".xlsx")
        )
    tables$comps_wide %>% 
        map_dfr(rbind, .id = "response") %>% 
        writexl::write_xlsx(
            paste0("export/models_comp_wide", Sys.Date(), ".xlsx")
        )
    
    # models pics text
    plots$text %>% 
        `[`(-1) %>% 
        map2(
            names(.), 
            ~ggsave(
                plot = .x, 
                filename = paste0("export/Fig.1_", .y, "_", Sys.Date(), ".png"), 
                height = 8, width = 11, dpi = 600)
        )
    
    # models pics suppl.
    plots$suppl %>% 
        `[`(-1) %>% 
        map2(
            names(.), 
            ~ggsave(
                plot = .x, 
                filename = paste0("export/Suppl.1_", .y, "_", Sys.Date(), ".png"), 
                height = 8, width = 11, dpi = 600)
        )
}

pdf("export/multipage_plots.pdf")
for(i in c("abundance", "abundance_log", "nsp", "nsp100", "shan")){
    res[[i]] %>% 
        select(-fits, -pred) %>% 
        mutate(aic = round(aic, 1), r2 = round(r2, 2)) %>% 
        tableGrob(rows = NULL) %>% 
        grid.arrange(plots[[i]], ., ncol = 1, heights = c(2,1))
}
dev.off()


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
