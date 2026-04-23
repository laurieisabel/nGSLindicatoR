# Fonction créée 16 sept. 2021 par Jordan Ouellette-Plante
#   modifée par jmc 

#' Fonction calculant les PUE/NUE ainsi que l'abondance/biomasse et les intervalles de confiance. Fournit les valeurs avec ou sans
#' le modèle multiplicatif.
#' 
#' Se base sur les scripts SAS retrouvés ici:
#'  - S:/SAS/Paces/Équivalent_Teleost/Program/pue.sas
#'  - S:/MissionNGSL/Analyse/ProgramSAS/pue_2.sas
#' Cette fonction est utilisée en parallèle à celle read_NGSL() qui fournit les intrants.

#' @param ans Vecteur des années disponibles ou à considérer.
#' @param set Dataframe des données de traits (i.e. set)
#' @param catch Dataframe des données de captures (i.e. catch)
#' @param stratum Dataframe des données de strates (i.e. stratum)
#' @param strates_sel Par défaut NULL. Fournir un vecteur de numéros de strates si on veut utiliser un sous-ensemble des strates. 
#' Par exemple, quelqu'un pourrait être intéressé à n'utiliser que les strates profondes du golfe.

#' @param compa.year.meth Méthode de calcul pour les années avec stations comparatives: 
#'                         Par défaut 1 (prends tous les traits comme s'ils étaient indépendants). 
#'                         2: fait la moyenne des stations comparatives avant de faire la moyenne de la strate
#'                         3: enlève les stations comparatives faites avec l'ancien navire 
#' @param ouv_hori_stan  Largeur d'ouverture du standard. Ex 16.94 teleost, 16.71 cabot
#' @param option_ouv_hori_stan options de calcul, option 1 est celle utilsée pour éval de stock 2022 soit celle à utliser avec les facteurs de conversion à Bourdages et al. 2007, option 2 est la "bonne" à utiliser avec les facteurs de conversion à Hugues
#' @keywords NGSl, catch, set, carbio

require(tidyverse)
options(dplyr.summarise.inform = FALSE)

get_pop_indices <- function(ans, set, catch, stratum, strates_sel = NULL, 
                            compa.year.meth = 3, ouv_hori_stan = 16.71,
                            option_ouv_hori_stan = 2){
  
  
  # 1. surface standard d'un trait effectué par le tandem navire--engin Teleost--Campelen ####
  a <- ouv_hori_stan  * 0.75 * 1852 / 10^6 # en km2
  # 16.94 m = ouverture horizontale
  # 0.75 MN = distance chalutée durant 15 minutes à 3 Noeuds/h (durée et vitesse visées avec le Teleost--Campelen)
  # 1852 = nbre de mètre dans un mille nautique
  
  # 2 tri des strates, au besoin ####
  if (is.null(strates_sel)){ # si rien de spécifié, on travaille avec les 56 strates normalement utilisées (celles ± échantillonnées au fil des années)
    strates_sel <- c(401,402,403,404,405,406,407,408,409,410,411,412,413,414,801,802,803,
                     804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,
                     821,822,823,824,827,828,829,830,831,832,833,835,836,837,838,839,840,
                     841,851,852,854,855)
  }
  stratum <- stratum %>% filter(no_str %in% strates_sel)
  
  # 3. Nbre d'unités chalutables par strates ####
  stratum <- stratum %>% mutate(nb_unit = surfkm91 / a)
  
  # 4. Sélection et préparation des données ####
  # données set: années voulus, traits aléatoires réussis
  set <- set %>% 
    filter(annee %in% ans, typtrait == 1, resultat %in% 1:2, no_str %in% strates_sel) %>% 
    left_join(., stratum, by = "no_str") %>% # ajout des superficies des strates / nb d'unités chalutables
    select(annee, nav, rel, trait, opano, no_str, dist_vit, ouv_hori, surfkm91, nb_unit) %>% # variables requises
    mutate(trait = paste0(nav, trait)) # pour distinger les traits d'une année comparative où il y aurait 2 navires
  
  # données catch: standardisation des captures pour un trait standard en équivalents NGCC Teleost--Campelen
  if(option_ouv_hori_stan == 1){
    catch <- catch %>% 
      mutate(trait = paste0(nav, trait)) %>% # pour distinger les traits d'une année comparative où il y aurait 2 navires
      left_join(., set %>% select(nav, rel, trait, dist_vit, ouv_hori), by = c("nav", "rel", "trait")) %>% 
      filter(!is.na(dist_vit)) %>% # si dist_vit est NA, la ligne provient d'un trait qu'on ne garde pas
      mutate(surf_ech = ouv_hori * (dist_vit * 1852) / 10^6, # en km2
             pds_capt_cor = pds_capt_cor_ca * a / surf_ech, # standardisation du trait
             nb_capt_cor = nb_capt_cor_ca * a / surf_ech) %>% # standardisation du trait
      group_by(nav, rel, trait) %>% 
      summarise_at(.vars = c("nb_capt_cor", "pds_capt_cor"), .funs = sum) %>% # combine les catégories, si présentes
      ungroup
  }

  if(option_ouv_hori_stan == 2){  # comme les facteurs de conversion à hugues tiennent compte de l'ouverture, on met l'ouverture standard pour tout le monde et on standardise juste avec distance
    catch <- catch %>% 
      mutate(trait = paste0(nav, trait)) %>% # pour distinger les traits d'une année comparative où il y aurait 2 navires
      left_join(., set %>% select(nav, rel, trait, dist_vit, ouv_hori), by = c("nav", "rel", "trait")) %>% 
      filter(!is.na(dist_vit)) %>% # si dist_vit est NA, la ligne provient d'un trait qu'on ne garde pas
      mutate(surf_ech = ouv_hori_stan * (dist_vit * 1852) / 10^6, # en km2
             pds_capt_cor = pds_capt_cor_ca * a / surf_ech, # standardisation du trait
             nb_capt_cor = nb_capt_cor_ca * a / surf_ech) %>% # standardisation du trait
      group_by(nav, rel, trait) %>% 
      summarise_at(.vars = c("nb_capt_cor", "pds_capt_cor"), .funs = sum) %>% # combine les catégories, si présentes
      ungroup
  }
  

  ### moyenner la capture des traits comparatifs
  if(compa.year.meth == 2){
    set$trait2 <- substr(set$trait,3,nchar(set$trait))# refaire variable trait pcq modifiée pour inclure numero de navire plus huat...
    set$annee_tra <- paste(set$annee, set$trait2, sep = "_")
    catch.comp <- merge(x=catch, y=set, all.x=T, all.y=F, by = c("nav", "rel", "trait"))
    catch.comp.agg <- aggregate(cbind(pds_capt_cor, nb_capt_cor) ~ annee_tra, data = catch.comp, mean) # moyenner pour stations comparatives
    #### amener info manquantes dans catch.comp.agg pour poursuivre
    catch.comp.agg2 <- merge(x = catch.comp.agg, 
                             y = set[!duplicated(set$annee_tra),], 
                             all.x=T, all.y = F, by = "annee_tra")
    catch2 <- catch.comp.agg2[,c('nav','rel','trait','nb_capt_cor','pds_capt_cor')]
    
    #### identifier et enlever les traits remplacés par la moyenne dans set, sinon capture sera de 0 à l'étape 5, ce qui fausse les résultats
    catch$nav_rel_tra <- paste(catch$nav, catch$rel, catch$trait, sep = "_")
    catch2$nav_rel_tra <- paste(catch2$nav, catch2$rel, catch2$trait, sep = "_")
    to.remove <- catch[!(catch$nav_rel_tra %in% catch2$nav_rel_tra),]$nav_rel_tra
    set$nav_rel_tra <- paste(set$nav, set$rel, set$trait, sep = "_")
    
    set <- set[!(set$nav_rel_tra %in% to.remove),]
    catch <- catch2
    
    rm(catch.comp, catch.comp.agg, catch.comp.agg2, to.remove, catch2)
  }
  
  
  # 5. Création de setcatch ####
  setcatch <- left_join(set, catch, by = c("nav", "rel", "trait")) %>% 
    mutate_at(.vars = c("nb_capt_cor", "pds_capt_cor"), 
              .funs = function(x) ifelse(is.na(x), 0, x)) # Les valeurs NA sont en fait des 0
  
  ## enlever les traits comparatifs si duplicata
  if(compa.year.meth == 3){
    ### identifier les traits comparatifs
    toto <- setcatch
    toto$trait2 <- substr(toto$trait,3,nchar(toto$trait))
    toto$annee_trait <- paste(toto$annee, toto$trait2, sep = "_")
    trait.comp <- toto[duplicated(toto$annee_trait),]$annee_trait
    # table(substr(trait.comp,1,4)) # 1990, 2004:2005,2021:2022
    
    ### identifier les traits à enlever: traits comparatifs "ancien" bateau
    toto$keep <- 1
    toto$keep <- ifelse(toto$annee == 1990 & toto$annee_trait %in% trait.comp & toto$nav == 31, 0, 
                        ifelse(toto$annee %in% c(2004:2005) & toto$annee_trait %in% trait.comp & toto$nav == 34, 0,
                               ifelse(toto$annee %in% 2022 & toto$annee_trait %in% trait.comp & toto$nav == 39, 0, toto$keep)))
    # sum(toto$keep== 0) # on enlève 249 stations
    col.to.keep <- c("annee","nav","rel","trait","opano","no_str","dist_vit","ouv_hori","surfkm91","nb_unit","nb_capt_cor","pds_capt_cor")
    setcatch <- toto[toto$keep == 1, col.to.keep]
    rm(toto, trait.comp)
  }
  
  # 6. Surface échantillonnée annuellement ####
  a_t0 <- set %>% 
    group_by(annee, no_str, surfkm91) %>% 
    count %>% # nbre de traits dans chaque cellule année-strate
    ungroup %>% 
    filter(n >= 1) %>% # ne conserve que les cellules année-strate où j'ai >= 2 traits
    # RAISON: ça me prend au moins 2 traits pour calculer une moyenne ou variance
    group_by(annee) %>% 
    summarise(a_t0 = sum(surfkm91)) %>% 
    ungroup
  
  # 7. Calculs sans le modèle multiplicatif ####
  # Calcul du nombre d'observations, de la moyenne et de la variance par année-strate
  tmp <- setcatch %>% 
    group_by(annee, no_str) %>% 
    mutate(n_traits = n_distinct(trait)) %>% 
    ungroup %>% 
    filter(n_traits >= 1) %>% # ne conserve que les cellules année-strate où j'ai >= 2 traits
    group_by(annee, no_str) %>% 
    summarise(n = n(),
              mean_n = mean(nb_capt_cor),
              mean_p = mean(pds_capt_cor),
              var_n = sd(nb_capt_cor),
              var_p = sd(pds_capt_cor),
              nb_unit = unique(nb_unit),
              surfkm91 = unique(surfkm91),
              total_n = mean_n * nb_unit, # nombre total de la strate
              total_p = mean_p * nb_unit) %>%  # biomasse totale de la strate
    ungroup
  
  # Calcul de la moyenne et de la variance pour la zone d'étude, selon un plan d'éch. aléatoire stratifié
  tmp <- left_join(tmp, a_t0, by = "annee") %>% 
    mutate(w = surfkm91 / a_t0, # poids de la strate p/r à l'aire d'étude
           wmeanp = w * mean_p,
           #wvarp = (w^2 * var_p) * ((nb_unit - n) / nb_unit) / n,
           wvarp = w * var_p,
           wmeann = w * mean_n,
           #wvarn = (w^2 * var_n) * ((nb_unit - n) / nb_unit) / n,
           wvarn = w * var_n,
           ahp = nb_unit * (nb_unit - n) / n,
           ahn = nb_unit * (nb_unit - n) / n,
           d1p = ahp * var_p,
           d1n = ahn * var_n,
           d2p = (ahp * var_p)^2 / (n - 1),
           d2n = (ahn * var_n)^2 / (n - 1)) %>% 
    group_by(annee) %>% # calcul à l'échelle de l'aire d'étude
    summarise(mean_p = sum(wmeanp),
              mean_n = sum(wmeann),
              var_p = sum(wvarp),
              var_n = sum(wvarn),
              d1p = sum(d1p),
              d1n = sum(d1n),
              d2p = sum(d2p),
              d2n = sum(d2n)) %>% 
    ungroup %>% 
    left_join(., a_t0, by = "annee") %>% 
    mutate(nb_unit_total = a_t0 / a,
           d_p = d1p^2 / d2p, # degrés de liberté, poids
           d_n = d1n^2 / d2n, # degrés de liberté, nbre
           t_p = qt(p = 0.975, df = d_p), # t de Student, poids
           t_n = qt(p = 0.975, df = d_n), # t de Student, nbre
           nb_mean = mean_n,
           pds_mean = mean_p,
           sd_p = var_p,
           surf_ech = a_t0,
           pds_icinf = mean_p - t_p * var_p^0.5, # I.C., borne inf., poids moyen
           pds_icsup = mean_p + t_p * var_p^0.5, # I.C., borne sup., poids moyen
           nb_icinf = mean_n - t_n * var_n^0.5, # I.C., borne inf., nbre moyen
           nb_icsup = mean_n + t_n * var_n^0.5, # I.C., borne sup., nbre moyen
           bimt_p = mean_p * nb_unit_total, # Biomasse totale de la zone (kg)
           bimt_n = mean_n * nb_unit_total, # Nombre totale de la zone
           bimt_p = pds_icinf * nb_unit_total, # I.C., borne inf., poids total
           bsmt_p = pds_icsup * nb_unit_total, # I.C., borne sup., poids total
           bimt_n = nb_icinf * nb_unit_total, # I.C., borne inf., nombre total
           bsmt_n = nb_icsup * nb_unit_total) %>% # I.C., borne sup., nombre total
    select(annee, surf_ech, nb_mean, nb_icinf, nb_icsup, pds_mean, sd_p, pds_icinf, pds_icsup)
  
  pue <- tmp
  
  # 10. Moyennes de la série
  
  ref <- pue %>%
    filter(annee %in% (min(ans):(max(ans)-1))) %>%  # la dernière année du relevé ne sert pas
    summarise(nb_ref_mean = mean(nb_mean),
              pds_ref_mean = mean(pds_mean),
              nb_ref_var = var(nb_mean),
              pds_ref_var = var(pds_mean)) %>%
    mutate(nb_ref_inf = nb_ref_mean - 0.5 * nb_ref_var^0.5,
           nb_ref_sup = nb_ref_mean + 0.5 * nb_ref_var^0.5,
           pds_ref_inf = pds_ref_mean - 0.5 * pds_ref_var^0.5,
           pds_ref_sup = pds_ref_mean + 0.5 * pds_ref_var^0.5)
  
  # 11. Configuration finale des données
  
  pue <- pue %>%
    mutate(nb_ref_mean = ref$nb_ref_mean,
           nb_ref_inf = ref$nb_ref_inf,
           nb_ref_sup = ref$nb_ref_sup,
           pds_ref_mean = ref$pds_ref_mean,
           pds_ref_inf = ref$pds_ref_inf,
           pds_ref_sup = ref$pds_ref_sup) %>%
    select(annee, surf_ech, nb_mean, nb_icinf, nb_icsup,
           pds_mean, sd_p, pds_icinf, pds_icsup,
           nb_ref_mean, nb_ref_inf, nb_ref_sup,
           pds_ref_mean, pds_ref_inf, pds_ref_sup) %>%
    # mutate_at(.vars = c("abd_icinf", "abd_icsup", "bio_icinf", "bio_icsup",
    #                     "nb_icinf", "nb_icsup", "bio_icinf", "bio_icsup",
    #                     "nb_icinf_0", "nb_icsup_0", "pds_icinf", "pds_icsup",
    #                     "pds_icinf_0", "pds_icsup_0"), 
    #           .funs = function(x) ifelse(x < 0, 0, x)) %>% # pour ne pas avoir des axes des y qui vont dans le -
    ## edit jmc: on veut que les intervalles de conf inf aillent sous zéro si c,est le cas, 
    ##   ça démontre que méthode paramétrique pas adéquate pour calculer intervalles de confiances dans certains cas. 
    mutate_at(.vars = c("surf_ech"),
              .funs = function(x) round(x, 1)) %>%
    mutate_at(.vars = c("nb_mean", "nb_icinf", "nb_icsup", "nb_ref_mean", "nb_ref_inf",
                        "nb_ref_sup", "pds_ref_mean", "pds_ref_inf", "pds_ref_sup"),
              .funs = function(x) round(x, 4)) %>%
    mutate_at(.vars = c("pds_mean", "sd_p", "pds_icinf", "pds_icsup"),
              .funs = function(x) round(x, 7))
  
  return(pue)    
}
