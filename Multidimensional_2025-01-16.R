# Data load ---------------------------------------------------------------
library(tidyverse)
theme_set(
    theme_bw() +
        theme(
            text = element_text(family = "serif", size = 15),
            legend.position = "bottom"
        )
)
L <- c("фон_2009", "фон_2014", "буф_2009", "буф_2014", "имп_2009", "имп_2014")
long <- readxl::read_excel("Data/Carabidae_25.01.2023.xlsx", 
                           sheet = "main_data") %>% 
    select(-duration, -traps) %>% 
    arrange(desc(site), taxa) %>% # -taxa
    mutate(tur = factor(tur), 
           site = fct_inorder(site),
           zone = case_when(zone == "fon" ~ "фоновая", 
                            zone == "bufer" ~ "буферная",
                            zone == "superimpact" ~ "суперимпактная", 
                            TRUE ~ "импактная"),
           zone = factor(zone, levels = c("фоновая", "буферная", "импактная", "суперимпактная")),
           km = str_extract(site, "[:digit:]{1,}"),
           km = as.numeric(km), 
           year = as.factor(year), 
           .after = "site")

# 1 = turs 1 and 2 are united
wide1 <- long %>% # div2 - turs are united
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    select(-no_insects)

dis <- wide1 %>% 
    mutate(
        zone = substr(zone, 1, 3),
        zone = fct_collapse(
        zone, `имп` = c("имп", "суп"))
        ) %>%
    select(-site) %>% 
    # group_by(year, zone, km, plot) %>% 
    # summarise_all(mean) %>% 
    # ungroup %>% 
    # # select_if(~ !is.numeric(.) || sum(.) != 0)
    unite("ID", zone, year, km, plot, sep = "_") %>% 
    # filter(rowSums(.[,2:ncol(.)]) > 0 
    #        # str_detect(ID, "супер", negate = TRUE)
    # ) %>%
    # filter(str_detect(ID, "супер", negate = TRUE))
    column_to_rownames("ID") %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
pc <- ape::pcoa(dis)
eig <- pc$values$Eigenvalues
eig <- round(eig/sum(eig)*100, 1)
pc <- pc$vectors %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    as_tibble()

# Ordination viz ----------------------------------------------------------
pc_data <- pc %>% 
    separate(ID, into = c("zone", "year", "site", "plot"), 
             sep = "_") %>% 
    mutate(#zone = str_replace_all(zone, "ая", "ый"),
           zone = factor(zone, levels = c("фон", "буф", "имп")))
           # zone = fct_relabel(zone, ~paste0(.x, " территория"))) 

pc_data %>% 
ggplot(aes(x = Axis.1, y = Axis.2, linetype = year,
               fill = zone, color = zone, shape = year)) +
    geom_point(color = "black", size = 2.5) +  
    stat_ellipse() + 
    scale_shape_manual(values = c(21, 22))+
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_fill_manual(values = c("green", "orange", "red")) +
    scale_color_manual(values = c("green", "orange", "red")) +
    labs(subtitle = "Ординация по динамической плотности\nТуры объединены", 
         x = paste0("Ось 1 (", eig[1], " %)"), 
         y = paste0("Ось 2 (", eig[2], " %)"), 
         # fill = "Год", color = "Год", shape = "Год"
         ) +
    theme(panel.grid = element_blank())
# ggsave("Plot5.svg", width = 18, height = 13, units = "cm")

# Distances ---------------------------------------------------------------
distances <- list(
    m1_raw = dis,
    m2_2axes = pc %>% 
        column_to_rownames("ID") %>%
        select(1:2) %>% 
        dist(),
    m3_allaxes = pc %>% 
        column_to_rownames("ID") %>% 
        dist()
) %>% 
    lapply(function(m){
        m <- as.matrix(m)
        m[upper.tri(m)] <- NA
        diag(m) <- NA
        m %>% 
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
    })
distances %>% 
    map(~mutate_if(.x, is.numeric, function(a){round(a, 2)})) %>% 
    writexl::write_xlsx("export/distances.xlsx")

# Permanova ---------------------------------------------------------------
permanova <- vegan::adonis2(
    dis ~ zone + year, data = wide1, 
    by = "terms", permutations = 9999)
permanova %>% 
    rownames_to_column("term") %>% 
    writexl::write_xlsx( "permanova.xlsx")

