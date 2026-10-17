# carga base original
RHC <- readRDS("RHC.rds")


# adm datos cáncer gastrico
RHC_CG <- RHC %>% 
  dplyr::filter(TOPOGRAFIA_CIE10 == "C16",
                EDAD > 29) %>% 
  dplyr::mutate(
    QUINQUENIO_EDAD  = cut(EDAD,
                           breaks = c(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, Inf),
                           labels = c( "30-34", "35-39", "40-44", "45-49","50-54", 
                                       "55-59", "60-64", "65-69", "70-74", "75-79", 
                                       "80 y más"),
                           right = FALSE))


# reconstrucción estadios == "99"
n_RHC_CG <- RHC_CG %>% 
  mutate(ESTADIO = as.character(ESTADIO)
  ) %>% 
  
  filter(
    
    # mantener solo casos con ESTADIO == 99
    ESTADIO        == "99",        
    
    # quitar filas con información desconocida en todas las variables TNM
    !(T_TUMOR        == "X" &        
        N_GANGLIO    == "X" & 
        M_METASTASIS == "X"),
    
    # quitar filas sin información en todas las variables TNM
    !(T_TUMOR        == "99" &        
        N_GANGLIO    == "99" & 
        M_METASTASIS == "99"),
    
  ) %>%              
  
  # reconstruir estadios a partir de M_METASTASIS == "1" -> hay metastasis
  mutate(n_ESTADIO = 
           if_else(M_METASTASIS %in% "1","IV", ESTADIO)
  ) %>%
  
  # mantener filas con n_ESTADIO  == "99" 
  filter(n_ESTADIO == "99") %>%
  
  # agrupar por combinaciones de TNM presente en el df
  group_by(T_TUMOR, N_GANGLIO, M_METASTASIS) %>%
  
  # número de casos por cada combinación existente
  summarise(Frecuencia = n(), .groups = "drop") %>%
  
  # ordenar combinaciones según frecuencia
  arrange(desc(Frecuencia), T_TUMOR, N_GANGLIO)