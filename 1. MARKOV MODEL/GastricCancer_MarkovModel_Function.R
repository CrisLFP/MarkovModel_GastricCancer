## Modelo de progredicón CG - FONIS

runs_gc_markov <- function(v_params, return_trace = FALSE) {
  
  with(as.list(v_params), {
  
  # create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle 
  
  ## creción multiple de 6 matrices (una por cada estrategia)
  m.M_no_prev <-  matrix(NA, ## llenar con valores NA
                                                                                         
                   ## row number n.t+1 = 83                                   
                   nrow = n.t + 1,  
                                                                                         
                   ## columns number (disease states) 
                   ncol = n.s, 
                                                                                         
                   ## cycles number
                   dimnames = list(paste("cycle", 0:n.t, sep = " "), 
                                   
                   ## disease states names
                   v.n))  
  
  

  ## select first row of each strategy matrix 
  m.M_no_prev[1, ] <-  c(
    
    ## population proportion at each state at first cycle
    ## Tabla 1 Ferrecio, et al 2009
    1 - 0.62,         # Normalhpn (sano, Hp-)
    0.62,             # Normalhpp (sano, Hp+)
    rep(0, n.s - 2)   # todos los demás estados en 0
  )   
  
## first row at matrix  
  m.M_no_prev[1, ]            <- 0                   ## todos los estados en 0
  m.M_no_prev[1, "Normalhpn"] <- 1 - 0.62            ## hp- en ciclo 1
  m.M_no_prev[1, "Normalhpp"] <- 0.62                ## hp+ en ciclo 1
  stopifnot(abs(sum(m.M_no_prev[1, ]) - 1) < 1e-12)  ## confirmar que en ciclo 1 todo sume 1
  
  
## ------- Array Transition Probabilities ------- ##
  
### Create the Array Transition Probabilities
  
  
  ##---------- Create Array Transition No Str ----------##
  
  a.P_noprev <- array(0,                                ## llenado incial con 0 en todo el array
                      
                      ## array dimmentions (82x39x39)
                      dim = c(n.s,                      ## 39 estados de salud
                              n.s,                      ## 39 estados de salud
                              n.t),                     ## 82 ciclos   
                      
                      ## labels
                      dimnames = list( from  = v.n,     ## nombres de estados
                                       to    = v.n,     ## nombres de estados
                                       cycle = 1:n.t))  ## número de ciclo
  
  
## ------- No Prevention strategy ------- ##
  
  
###  fill in the transition probability array no prevention
## EL ARRAY DE TRANSICION TIENE 70 CICLOS, PERO AHORA, MI MAXIMO DE CICLOS ES 82
  
  ### From Normal helicobacter (-)
  ## filas de "Normalhpn",  "Normalhpn", esta llenado con el objeto indicado
  a.P_noprev["Normalhpn",  "Normalhpn",]            <- 1 - p.NnD - p.NnNp - p.NnGn   ## quedarse en NHp-
  a.P_noprev["Normalhpn", "Normalhpp", ]            <- p.NnNp                        ## pasar de Nhp- a Nhp+
  a.P_noprev["Normalhpn", "Gastritishpn", ]         <- p.NnGn                        ## pasar  de Nhp- a GHp-
  a.P_noprev["Normalhpn", "Dead", ]                 <- p.NnD                         ## pasar de Hp- a muerte
  
  ### From Normal helicobacter (+) 
  a.P_noprev["Normalhpp", "Normalhpp", ]            <- 1 - p.NpGp - p.NpNn - p.NpD   ## quedarse en NHp+
  a.P_noprev["Normalhpp", "Gastritishpp", ]         <- p.NpGp                        ## pasar de NHp+ a GHp+
  a.P_noprev["Normalhpp", "Normalhpn", ]            <- p.NpNn                        ## pasar e NHp+ a Nhp-
  a.P_noprev["Normalhpp", "Dead", ]                 <- p.NpD                         ## pasar de NHp+ a muerte
  
  ### From Gastritis Helicobacter (-)
  a.P_noprev["Gastritishpn", "Gastritishpn", ]      <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD  ## quedarse en GHp-
  a.P_noprev["Gastritishpn", "Atrophyhpn", ]        <- p.GnAn                                ## pasar de GHp- a Atr Hp-
  a.P_noprev["Gastritishpn", "Gastritishpp", ]      <- p.GnGp                                ## pasar de GHp- a GHp+
  a.P_noprev["Gastritishpn", "Normalhpn", ]         <- p.GnNn                                ## pasar de GHp- a NHp-
  a.P_noprev["Gastritishpn", "Dead", ]              <- p.GnD                                 ## pasat fr GHp- a muerte
  
  ### From Gastritis Helicobacter (+)
  a.P_noprev["Gastritishpp", "Gastritishpp", ]      <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD  ## quedarse en GHp+
  a.P_noprev["Gastritishpp", "Atrophyhpp", ]        <- p.GpAp                                ## pasar de GHp+ a Atr Hp+
  a.P_noprev["Gastritishpp", "Gastritishpn", ]      <- p.GpGn                                ## pasar de GHp+ a GHp-
  a.P_noprev["Gastritishpp", "Normalhpp", ]         <- p.GpNp                                ## pasar de GHp+ a NHp+
  a.P_noprev["Gastritishpp", "Dead", ]              <- p.GpD                                 ## pasar de GHp+ a muerte
  
  ### From Atrophy Helicobacter (-) 
  a.P_noprev["Atrophyhpn", "Atrophyhpn", ]          <- 1 - p.AI - p.AnAp - p.AnGn - p.AnD    ## quedarse en Atr Hp-
  a.P_noprev["Atrophyhpn", "Intestinalhpn", ]       <- p.AI                                  ## pasar de Atr Hp- a IHp-
  a.P_noprev["Atrophyhpn", "Atrophyhpp", ]          <- p.AnAp                                ## pasar de Atr Hp- a Atr Hp+
  a.P_noprev["Atrophyhpn", "Gastritishpn", ]        <- p.AnGn                                ## pasar de Atr Hp- a GHp-
  a.P_noprev["Atrophyhpn", "Dead", ]                <- p.AnD                                 ## pasar de Atr Hp- a muerte
  
  ### From Atrophy Helicobacter (+) 
  a.P_noprev["Atrophyhpp", "Atrophyhpp", ]          <- 1 - p.AI - p.ApAn - p.ApGp - p.ApD    ## quedarse en Atr Hp+
  a.P_noprev["Atrophyhpp", "Intestinalhpp", ]       <- p.AI                                  ## pasar de Atr Hp+ a IHp+
  a.P_noprev["Atrophyhpp", "Atrophyhpn", ]          <- p.ApAn                                ## pasar de Atr Hp+ a Atr Hp-
  a.P_noprev["Atrophyhpp", "Gastritishpp", ]        <- p.ApGp                                ## pasar de Atr Hp+ a GHp+
  a.P_noprev["Atrophyhpp", "Dead", ]                <- p.ApD                                 ## pasar de Atr Hp+ a muerte
  
  ### From Intestinal Helicobacter (-)
  a.P_noprev["Intestinalhpn", "Intestinalhpn", ]    <- 1 - p.ID - p.InIp - p.InAn - p.InD    ## quedarse en I Hp-
  a.P_noprev["Intestinalhpn", "Dysplasiahpn", ]     <- p.ID                                  ## pasar de I Hp- a Dis Hp-
  a.P_noprev["Intestinalhpn", "Intestinalhpp", ]    <- p.InIp                                ## pasar de I Hp- a I Hp+
  a.P_noprev["Intestinalhpn", "Atrophyhpn", ]       <- p.InAn                                ## pasar de I Hp- a Atr Hp-
  a.P_noprev["Intestinalhpn", "Dead", ]             <- p.InD                                 ## pasar de I Hp- a muerte
  
  ### From Intestinal Helicobacter (+)
  a.P_noprev["Intestinalhpp", "Intestinalhpp", ]    <- 1 - p.ID - p.IpIn - p.IpAp - p.IpD    ## quedarse en I Hp+
  a.P_noprev["Intestinalhpp", "Dysplasiahpp", ]     <- p.ID                                  ## pasar de I Hp+ a Dys Hp+
  a.P_noprev["Intestinalhpp", "Intestinalhpn", ]    <- p.IpIn                                ## pasar de I Hp+ a I Hp+
  a.P_noprev["Intestinalhpp", "Atrophyhpp", ]       <- p.IpAp                                ## pasar de I Hp+ Atr Hp-
  a.P_noprev["Intestinalhpp", "Dead", ]             <- p.IpD                                 ## pasar de I Hp+ a muerte
  
  ### From Dysplasia Helicobacter (-)
  a.P_noprev["Dysplasiahpn", "Dysplasiahpn", ]      <- 1 - p.DnDp - p.DnIn - p.D1p - p.DnD   ## quedarse en Dys Hp- 
  a.P_noprev["Dysplasiahpn", "Dysplasiahpp", ]      <- p.DnDp                                ## pasar de Dys Hp- a Dys Hp+
  a.P_noprev["Dysplasiahpn", "Intestinalhpn", ]     <- p.DnIn                                ## pasar de Dys Hp- a I Hp-
  a.P_noprev["Dysplasiahpn", "1preclinical", ]      <- p.D1p                                 ## pasar de Dys Hp- a preclinical 1
  a.P_noprev["Dysplasiahpn", "Dead", ]              <- p.DnD                                 ## pasar de Dys Hp- a muerte
  
  ### From Dysplasia Helicobacter (+)
  a.P_noprev["Dysplasiahpp", "Dysplasiahpp", ]      <- 1 - p.DpDn - p.DpIp - p.D1p - p.DpD   ## quedarse en Dys Hp+
  a.P_noprev["Dysplasiahpp", "Dysplasiahpn", ]      <- p.DpDn                                ## pasar de Dys Hp+ a Dys Hp-
  a.P_noprev["Dysplasiahpp", "Intestinalhpp", ]     <- p.DpIp                                ## pasar de Dys Hp+ a I Hp +
  a.P_noprev["Dysplasiahpp", "1preclinical", ]      <- p.D1p                                 ## pasar de Dys Hp+ a preclinical 1
  a.P_noprev["Dysplasiahpp", "Dead", ]              <- p.DpD                                 ## pasar de Dys Hp+ a muerte
  
  ### From 1 preclinical 	
  a.P_noprev["1preclinical", "1preclinical",]       <- 1 - p.1p1ca - p.1p2p - p.1pD	        ## quedarse en preclinico 1           
  a.P_noprev["1preclinical", "1clinicala",]         <- p.1p1ca	                            ## pasar de preclinico 1 clinico 1
  a.P_noprev["1preclinical", "2preclinical",]       <- p.1p2p	                              ## pasar de preclinico 1 a preclinico 2
  a.P_noprev["1preclinical", "Dead", ]              <- p.1pD	                              ## pasar de preclinico 1 muerte
  
  ### From 2 preclinical	
  a.P_noprev["2preclinical", "2preclinical",]       <- 1 - p.2p2ca - p.2p3p - p.2pD	        ## quedarse en preclinico 2
  a.P_noprev["2preclinical", "2clinicala", ]        <- p.2p2ca	                            ## pasar de preclinico 2 a clinico 2
  a.P_noprev["2preclinical", "3preclinical", ]      <- p.2p3p	                              ## pasar de preclinico 2 a preclinico 3
  a.P_noprev["2preclinical", "Dead", ]              <- p.2pD	                              ## pasar de preclinico a muerte
  
  ### From 3 preclinical	
  a.P_noprev["3preclinical", "3preclinical", ]      <- 1 - p.3p3ca - p.3p4p - p.3pD	        ## quedarse en preclinico 3         
  a.P_noprev["3preclinical", "3clinicala", ]        <- p.3p3ca	                            ## pasar de preclinico 3 clinico 3
  a.P_noprev["3preclinical", "4preclinical", ]      <- p.3p4p	                              ## pasar de preclinico 3 a clinico 4
  a.P_noprev["3preclinical", "Dead", ]              <- p.3pD	                              ## pasar de preclinico 3 a muerte
  
  ### From 4 preclinical	
  a.P_noprev["4preclinical", "4preclinical", ]      <- 1 - p.4p4ca - p.4pD	                ## quedarse en preclinico 4                          
  a.P_noprev["4preclinical", "4clinicala", ]        <- p.4p4ca	                            ## pasar de preclinico 4 a clinco 4
  a.P_noprev["4preclinical", "Dead", ]              <- p.4pD	                              ## pasar de preclinico 4 a muerte
  
  ### From 1 clinical a	
  a.P_noprev["1clinicala", "1clinicala", ]          <- 0	                                  ## quedarse en clinico 1a
  a.P_noprev["1clinicala", "1clinicalb", ]          <- 1 - p.1caD	                          ## pasar de clinico 1a a 1b
  a.P_noprev["1clinicala", "Dead", ]                <- p.1caD	                              ## pasar de clinico 1a a muerte
  
  ### From 2 clinical a	
  a.P_noprev["2clinicala", "2clinicala", ]          <- 0	                                  ## quedarse en clinico 2a
  a.P_noprev["2clinicala", "2clinicalb", ]          <- 1 - p.2caD	                          ## pasar de clinico 2a a 2b
  a.P_noprev["2clinicala", "Dead", ]                <- p.2caD	                              ## pasar de clinico 2a a muerte
  
  ### From 3 clinical a	
  a.P_noprev["3clinicala", "3clinicala", ]          <- 0 	                                  ## quedarse en clinico 3a
  a.P_noprev["3clinicala", "3clinicalb", ]          <- 1 - p.3caD		                        ## pasar de clinico 3a a 3b
  a.P_noprev["3clinicala", "Dead", ]                <- p.3caD	                              ## pasar de clinico 3a a muerte
  
  ### From 4 clinical a 	
  a.P_noprev["4clinicala", "4clinicala", ]          <- 0	                                  ## quedarse en clinico 4a
  a.P_noprev["4clinicala", "4clinicalb", ]          <- 1 - p.4caD  	                        ## pasar de clinico 4a a 4b
  a.P_noprev["4clinicala", "Dead", ]                <- p.4caD	                              ## pasar de clinico 4a a muerte
  
  ### From 1 clinical b	
  a.P_noprev["1clinicalb", "1clinicalb", ]          <- 0	                                  ## quedarse en clinico 1b
  a.P_noprev["1clinicalb", "1clinicalc", ]          <- 1 - p.1cbD	                          ## pasar de clinico 1b a 1c
  a.P_noprev["1clinicalb", "Dead", ]                <- p.1cbD	                              ## pasar de clinico 1b a muerte
  
  ### From 2 clinical b	
  a.P_noprev["2clinicalb", "2clinicalb", ]          <- 0	                                  ## quedarse en clinico 2b
  a.P_noprev["2clinicalb", "2clinicalc", ]          <- 1 - p.2cbD	                          ## pasar de clinico 2b a 2c
  a.P_noprev["2clinicalb", "Dead", ]                <- p.2cbD	                              ## pasar de clinico 2b a muerte
  
  ### From 3 clinical b	
  a.P_noprev["3clinicalb", "3clinicalb", ]          <- 0 	                                  ## quedarse en clinico 3b
  a.P_noprev["3clinicalb", "3clinicalc", ]          <- 1 - p.3cbD	                          ## pasar de clinico 3b a 3c
  a.P_noprev["3clinicalb", "Dead", ]                <- p.3cbD	                              ## pasar de clinico 3b a muerte
  
  ### From 4 clinical b	
  a.P_noprev["4clinicalb", "4clinicalb", ]          <- 0 	                                  ## quedarse en clinico 4b
  a.P_noprev["4clinicalb", "4clinicalc", ]          <- 1 - p.4cbD	                          ## pasar de clinico 4b a 3c
  a.P_noprev["4clinicalb", "Dead", ]                <- p.4cbD	                              ## pasar de clinico 4b a muerte
  
  ### From 1 clinical c	
  a.P_noprev["1clinicalc", "1clinicalc", ]          <- 0                                    ## qudarse en clinico 1c
  a.P_noprev["1clinicalc", "1clinicald", ]          <- 1 - p.1ccD	                          ## pasar de clinico 1c a 1d
  a.P_noprev["1clinicalc", "Dead", ]                <- p.1ccD	                              ## pasar de clinico 1c a muerte
  
  ### From 1 clinical d	
  a.P_noprev["1clinicald", "1clinicald", ]          <- 0                                   ## quedarse en clinico 1d
  a.P_noprev["1clinicald", "1clinicale", ]          <- 1 - p.1cdD	                         ## pasar de clinico 1d a 1e
  a.P_noprev["1clinicald", "Dead", ]                <- p.1cdD	                             ## pasar de clinico 1d a muerte
  
  ### From 1 clinical e	
  a.P_noprev["1clinicale", "1clinicale", ]          <- 0                                   ## quedarse en clinico 1e
  a.P_noprev["1clinicale", "1clinicalf", ]          <- 1 - p.1ceD	                         ## pasar de clinico 1e a 1f
  a.P_noprev["1clinicale", "Dead", ]                <- p.1ceD	                             ## pasar de clinico 1e a muerte
  
  ### From 1 clinical f	
  a.P_noprev["1clinicalf", "1clinicalf", ]          <- 1 - p.1cfNn - p.1cfD                ##
  a.P_noprev["1clinicalf", "Normalhpp", ]           <- p.1cfNn	                                       
  a.P_noprev["1clinicalf", "Dead", ]                <- p.1cfD	
  
  ### From 2 clinical c	
  a.P_noprev["2clinicalc", "2clinicalc", ]          <- 0
  a.P_noprev["2clinicalc", "2clinicald", ]          <- 1 - p.2ccD	
  a.P_noprev["2clinicalc", "Dead", ]                <- p.2ccD	
  
  ### From 2 clinical d	
  a.P_noprev["2clinicald", "2clinicald", ]          <- 0
  a.P_noprev["2clinicald", "2clinicale", ]          <- 1 - p.2cdD	
  a.P_noprev["2clinicald", "Dead", ]                <- p.2cdD	
  
  ### From 2 clinical e	
  a.P_noprev["2clinicale", "2clinicale", ]          <- 0
  a.P_noprev["2clinicale", "2clinicalf", ]          <- 1 - p.2ceD	
  a.P_noprev["2clinicale", "Dead", ]                <- p.2ceD	
  
  ### From 2 clinical f	
  a.P_noprev["2clinicalf", "2clinicalf", ]          <- 1 - p.2cfNn - p.2cfD                            
  a.P_noprev["2clinicalf", "Normalhpp", ]           <- p.2cfNn                                         
  a.P_noprev["2clinicalf", "Dead", ]                <- p.2cfD	
  
  ### From 3 clinical c	
  a.P_noprev["3clinicalc", "3clinicalc", ]          <- 0 	
  a.P_noprev["3clinicalc", "3clinicald", ]          <- 1 - p.3ccD
  a.P_noprev["3clinicalc", "Dead", ]                <- p.3ccD	
  
  ### From 3 clinical d	
  a.P_noprev["3clinicald", "3clinicald", ]          <- 0
  a.P_noprev["3clinicald", "3clinicale", ]          <- 1 - p.3cdD	
  a.P_noprev["3clinicald", "Dead", ]                <- p.3cdD	
  
  ### From 3 clinical e	
  a.P_noprev["3clinicale", "3clinicale", ]          <- 0
  a.P_noprev["3clinicale", "3clinicalf", ]          <- 1 - p.3ceD	
  a.P_noprev["3clinicale", "Dead", ]                <- p.3ceD	
  
  ### From 3 clinical f	
  a.P_noprev["3clinicalf", "3clinicalf", ]          <- 1 - p.3cfNn - p.3cfD                             
  a.P_noprev["3clinicalf", "Normalhpp", ]           <- p.3cfNn	                                        
  a.P_noprev["3clinicalf", "Dead", ]                <- p.3cfD	
  
  ### From 4 clinical c	
  a.P_noprev["4clinicalc", "4clinicalc", ]          <- 0 	
  a.P_noprev["4clinicalc", "4clinicald", ]          <- 1 - p.4ccD	
  a.P_noprev["4clinicalc", "Dead", ]                <- p.4ccD	
  
  ### From 4 clinical d	
  a.P_noprev["4clinicald", "4clinicald", ]          <- 0
  a.P_noprev["4clinicald", "4clinicale", ]          <- 1 - p.4cdD	
  a.P_noprev["4clinicald", "Dead", ]                <- p.4cdD	
  
  ### From 4 clinical e	
  a.P_noprev["4clinicale", "4clinicale", ]          <- 0
  a.P_noprev["4clinicale", "4clinicalf", ]          <- 1 - p.4ceD	
  a.P_noprev["4clinicale", "Dead", ]                <- p.4ceD	
  
  ### From 4 clinical f	
  a.P_noprev["4clinicalf", "4clinicalf", ]          <- 1 - p.4cfNn - p.4cfD                             
  a.P_noprev["4clinicalf", "Normalhpp", ]           <- p.4cfNn	                                        
  a.P_noprev["4clinicalf", "Dead", ]                <- p.4cfD	
  
  ### From Dead
  a.P_noprev["Dead", "Dead", ]                      <- 1 
  

##---------- Checks ----------## 
  
  rowSums(a.P_noprev) 
  
  ## validación por ciclo
  validacion_filas <- apply(a.P_noprev, MARGIN = 3, FUN = rowSums)
  
  any(a.P_noprev < 0)
  
  
### Identify negative probabilities 
  
  ## coordenates
  neg_probs_a.P_noprev <- which(a.P_noprev < 0,    ## condition
                                arr.ind = T)       ## cells
  
  ## df
  df_neg_probs_a.P_noprev <- data.frame(
    indx_inc_state = neg_probs_a.P_noprev[, "from"],
    indx_end_state = neg_probs_a.P_noprev[, "to"],
    cycle          = neg_probs_a.P_noprev[,"cycle"],
    value          = a.P_noprev[neg_probs_a.P_noprev]
  )
  
  ## state names
  df_neg_probs_a.P_noprev$inc_state <- v.n[df_neg_probs_a.P_noprev$indx_inc_state]
  df_neg_probs_a.P_noprev$end_state <- v.n[df_neg_probs_a.P_noprev$indx_end_state]
  
  ## order df
  df_neg_probs_a.P_noprev <- df_neg_probs_a.P_noprev[, c("cycle", "inc_state", "end_state", "value")]
  nrow(df_neg_probs_a.P_noprev)

  
  ## apply competitive_risk_function whe you have negative probabilites in the transition array
  a.P_noprev <- competitive_risk_correction(
    transition_array = a.P_noprev,
    df_neg_probs = df_neg_probs_a.P_noprev
  )

  ## check the re-scaled transition array
  rowSums(a.P_noprev) 
  validacion_filas <- apply(a.P_noprev, MARGIN = 3, FUN = rowSums)
  any(a.P_noprev < 0)

  
##---------- Markov Model ----------##
  
  
  ## cadena markov sin intervención 
  ## recorrer 82 ciclos, 1 ciclo equivale a 1 año
  for (t in 1:n.t){                                                   ## loop through the number of cycles. Desde t=1 hasta t = n.t
    
    ## en el ciclo t+1 se recibe lo que ocurrió en el ciclo t
    ## mutiplicar la distribuci´pon de la cochorte en el ciclo t (vectro de 1x39), por el array de transicion de ese ciclo
    m.M_no_prev[t + 1, ] <- m.M_no_prev[t, ] %*% a.P_noprev [ , , t]  ## estimate the Markov trace for cycle the next cycle (t + 1)
  }
  
  
  ## output of interest
  
  ## 1.def columns index
  ## grep() se devuelve la posicion en v.n de las columnas con los nombres
  
  ## normal stages
  idx_Nn <- grep("Normalhpn", v.n)       ## posicion normal hpn 
  idx_Np <- grep("Normalhpp", v.n)       ## posicion normal hpp
  
  ## pre cancer stages hp+
  idx_Gp <- grep("Gastritishpp",  v.n)  ## posicion gastritis hpp
  idx_Ap <- grep("Atrophyhpp",    v.n)  ## posicion atrofia hpp
  idx_Ip <- grep("Intestinalhpp", v.n)  ## posicion metapalsia hpp
  idx_Dp <- grep("Dysplasiahpp",  v.n)  ## posicion displasia hpp
  
  ## pre cancer stages hp-
  idx_Gn <- grep("Gastritishpn",  v.n)  ## posicion gastritis hpn
  idx_An <- grep("Atrophyhpn",    v.n)  ## posicion atrofia hpn
  idx_In <- grep("Intestinalhpn", v.n)  ## posicion metapalsia hpn
  idx_Dn <- grep("Dysplasiahpn",  v.n)  ## posicion displasia hpn
  
  ## cancer
  idx_clinical_all <- grep("clinical", v.n)       ## posicion  estadios neoplásicos (pre y clinicos)                              
  idx_preclinical  <- grep("preclinical", v.n)    ## posicion estadios preclinicos                              
  idx_cancer       <- setdiff(idx_clinical_all, idx_preclinical) ## posicion estadios cancer clinico 
  idx_dead         <- grep("Dead", v.n)                          ## posicion muerte
  idx_risk         <- setdiff(1:n.s, c(idx_cancer, idx_dead))    ## posicion estadios en riesgo de cancer
  
  ## cancer stages
  idx_stg1 <- grep("1clinical", v.n)  ## posicion clinico 1
  idx_stg2 <- grep("2clinical", v.n)  ## posicion clinico 2
  idx_stg3 <- grep("3clinical", v.n)  ## posicion clinico 3
  idx_stg4 <- grep("4clinical", v.n)  ## posicion clincio 4
  
  ## prevalencias marginales
  ## funcion para sumar la ocupacion de que tiene cada estado del denominador,
  ## en todas las filas del arg "rows"; y la divide por la suma de los estados del denominador en las mismas filas
  ## se obtiene una sola prevalencia, por algun grupo de filas (tramo etario)
  prev_marginal <- function(mat,
                            idx_estado, ## indices del numerador
                            idx_denom,  ## indices del denominador
                            rows) {     ## filas
    
    sum(mat[rows, idx_estado]) / sum(mat[rows, idx_denom])
  }
  
  ## prev atrophy
  idx_alive <- setdiff(1:n.s, idx_dead)                 ## todos los estadios vivos
  idx_hpp   <- intersect(grep("hpp", v.n), idx_alive)   ## todos los estadios con Hp+ vivos (quedan excluidos los pre y clinicos)
  idx_hpn   <- intersect(grep("hpn", v.n), idx_alive)   ## todos los estadios con Hp- vivos (quedan excluidos los pre y clinicos)
  
  ## age range fot evaluate nacional prevalence
  rows_prev <- 13:83                             ## filas 13 a 83 de la traza (30-100 años)
  
  
  ## pob 2022 - Pob P Zitko
  ## pop_ine_tramo <- c(6153716, 2695688, 2488314, 2222376, 1519974, 450362)  
  ## pop_ine_total <- 19828563 
  
  ## pob segi para estandarizar las tasas del modelo.
  ## pesos pob munial segi, agregados en tramos
  segi_w     <- c(40000, 28000, 12000, 9000, 7000, 3000, 1000)  
  segi_total <- sum(segi_w)
  
  ## 2. empty vector for n.t cycle
  ## creacion de vecroes con 82 ceros para guardar la proporcion de cada estadio por ciclo
  v_new_cases       <- numeric(n.t)
  v_cancer_deaths   <- numeric(n.t)
  v_pop_risk_cycle  <- numeric(n.t)
  v_pop_alive_cycle <- numeric(n.t)
  
  v_prop_stg1 <- numeric(n.t)
  v_prop_stg2 <- numeric(n.t)
  v_prop_stg3 <- numeric(n.t)
  v_prop_stg4 <- numeric(n.t)
  
  
  ## 3.1. calculates cycle by cycle
  for (t in 1:n.t) {          ## recorreo los ciclos 1:82, cada cilo t usa la fila t
    
    
    ## Clinical stages proportions:
    
      ## Denominador para estadios: Suma de todos los cánceres clínicos en el ciclo t 
      ## pop con cancer en ciclo t
      pop_cancer_clin_t <- sum(m.M_no_prev[t, idx_cancer])
      
      # Evitar división por cero si aún no hay casos de cáncer (ejemplo, edades tempranas)
      if(pop_cancer_clin_t > 0){
        
        ## prop de pacientes clincos que están en cada estadio en el ciclo t
        ## de todos los que tienen cancer, cuantos estan en cada estadios en cada ciclo
        v_prop_stg1[t] <- sum(m.M_no_prev[t, idx_stg1]) / pop_cancer_clin_t
        v_prop_stg2[t] <- sum(m.M_no_prev[t, idx_stg2]) / pop_cancer_clin_t
        v_prop_stg3[t] <- sum(m.M_no_prev[t, idx_stg3]) / pop_cancer_clin_t
        v_prop_stg4[t] <- sum(m.M_no_prev[t, idx_stg4]) / pop_cancer_clin_t
      } else {
        
        # sino hay casos, dejar las proporciones en 0
        v_prop_stg1[t] <- v_prop_stg2[t] <- v_prop_stg3[t] <- v_prop_stg4[t] <- 0
      }
    
    # population at each cycle (denominators)
    pop_risk_t           <- m.M_no_prev[t, idx_risk]     ## vector de poblacion sin cancer cada ciclo
    pop_alive_t          <- m.M_no_prev[t, -idx_dead]    ## vector de poblacion viva cada ciclo
    pop_cancer_t         <- m.M_no_prev[t, idx_cancer]   ## vector poblacion con cancer cada ciclo
    
    v_pop_risk_cycle[t]  <- sum(pop_risk_t)              ## total pob sin cancer en cada ciclo
    v_pop_alive_cycle[t] <- sum(pop_alive_t)             ## total pob viva en caca ciclo
    
    # incidence (risk pop * cancer prob)
    trans_prob_inc       <- a.P_noprev[idx_risk, idx_cancer, t]   ## matriz de probs del ciclo t, de estados en riesgo(preclinical) hacia cancer(clinical) 
    v_new_cases[t]       <- sum(pop_risk_t %*% trans_prob_inc)    ## casos nuevos a estadios clinicos en cada ciclo: pob en riesgo * prob de transicion a cancer
    
    # mortality (Pop with cáncer * death prob)
    prob_death_cancer    <- a.P_noprev[idx_cancer, idx_dead, t]      ## vector con las probs de morir en cada estadio clínico en cada ciclo
    v_cancer_deaths[t]   <- sum(pop_cancer_t %*% prob_death_cancer)  ## muertes nuevas de cancer en cada ciclo: pob con ancer * prob de morir
    
  }
  
    ## 3.2 prevalencias pre-cancer en ciclos específicos
    
    burn_in <- 12                             ## ciclos previos de loa 18 a 30 años
    ciclos_eval <- c(6, 21, 32, 46) + burn_in ## a partir del paper de Bellolio 2019.
    
    ## tramos de lesiones Bellolio
    l_filas_tramos_NG <- list(
      "19-44" = 2:27,
      "45-56" = 28:39,
      "57-66" = 40:49,
      "67-83" = 50:66
    )
    
    ## Calcular prev en esos puntos respecto a población viva
    ## Funcion que para cada fila del arg = rows, divide la suma de los estados idx, por la
    ## suma de de los vivos. Se ib un vector con la prevalencia por edad
    calc_prev <- function(mat, idx, rows) {
      
      num   <- rowSums(mat[rows, idx, drop = F])
      denom <- rowSums(mat[rows, -idx_dead])
      
      prevalencias <- num / denom
      return(as.numeric(prevalencias))
    }
    
    
    ## 3.3 Tasa de mortalidad por ciclo
    v_tasa_mort_ciclo <- (v_cancer_deaths / v_pop_alive_cycle) * 100000
    
    ## 3.4 tasa de incidencia por ciclo
    v_tasa_inc_ciclo <- (v_new_cases / v_pop_alive_cycle) * 100000
    
    ## 3.5 incidencia de pob en riesgo por ciclo
    v_tasa_inc_ciclo_riesgo <- (v_new_cases / v_pop_risk_cycle)  * 100000
    
    
    
  ## 4. Adjust for year 100 (cycle n.t + 1)
  ## Add 0 at end for align  83 rows (ages 18 to 100) -> recordar que la edad no transita a nada
  v_new_cases       <- c(v_new_cases, 0)
  v_cancer_deaths   <- c(v_cancer_deaths, 0)
  v_pop_risk_cycle  <- c(v_pop_risk_cycle, 0)
  v_pop_alive_cycle <- c(v_pop_alive_cycle, 0)
  
  ## 5. age group 
  ## Vector de 83 etiquetas: Ages 80-100: 1 final group of 21 years
  v_agrupacion <-c (
    rep(NA_real_, 2),    # años 18 y 19
    rep(1,        20),   # 20-39
    rep(2,        10),   # 40-49
    rep(3,        10),   # 50-59
    rep(4,        10),   # 60-69
    rep(5,        10),   # 70-79
    rep(6,        20),   # 80-99
    NA_real_             # edad 100      
  )
  
  ## stop contra desalineacion
  stopifnot(length(v_agrupacion) == length(v_new_cases))
  
  ## sum events and populations for each age group 
  ## sumar los casos nuevos por cada tramo. NA quedan excluidos
  v_inc_tramo       <- tapply(v_new_cases,       v_agrupacion, sum) 
  v_mort_tramo      <- tapply(v_cancer_deaths,   v_agrupacion, sum)
  v_pop_risk_tramo  <- tapply(v_pop_risk_cycle,  v_agrupacion, sum)
  v_pop_alive_tramo <- tapply(v_pop_alive_cycle, v_agrupacion, sum)
  
  
  ## 6. Calculate final rates per 100.000 habs 
  ## estas tasas son las análagas a las tasas brutas del modelo. no estan estandartizadas.
  v_tasa_mort_modelada <- (v_mort_tramo / v_pop_alive_tramo) * 100000   ## death respect alive pop at each cycle
  v_tasa_inc_modelada  <- (v_inc_tramo / v_pop_alive_tramo) * 100000     ## TARGET: cases respect alive pop 
  
  v_tasa_inc_riesgo    <- (v_inc_tramo / v_pop_risk_tramo)   * 100000  ## descriptive: cases respect risk pop at each cycle

  
  ## quite name if exists at this point
  ## quita los nombres arrastrados del tapply.
  v_tasa_inc_modelada  <- unname(v_tasa_inc_modelada)
  v_tasa_mort_modelada <- unname(v_tasa_mort_modelada)
  v_tasa_inc_riesgo    <- unname(v_tasa_inc_riesgo)
  
  
  # stopifnot(length(v_tasa_inc_modelada) == length(pop_ine_tramo))
  # tasa_inc_nacional <- sum(v_tasa_inc_modelada * pop_ine_tramo) / pop_ine_total
  stopifnot(length(v_tasa_inc_modelada) == length(segi_w) - 1L)   # 6 tramos del modelo (20-99)
  
  tasa_inc_nacional <- sum(c(0,                                   ## numerador 0 a menores de 20
                             v_tasa_inc_modelada) * segi_w) / segi_total
  
  
  ## ASR Segi de la tasa descriptiva (población en riesgo)
  tasa_inc_nacional_riesgo <- sum(c(0,                            ## numerador 0 a menores de 20
                                    v_tasa_inc_riesgo) * segi_w) / segi_total
  
  
  ## 7. proporciones por estadio agrgadas 
  
  ## 7.1 Descriptivo: suma total de pacientes en cada estadio a lo largo de todos los ciclos
  pop_total_stg1 <- sum(m.M_no_prev[, idx_stg1])
  pop_total_stg2 <- sum(m.M_no_prev[, idx_stg2])
  pop_total_stg3 <- sum(m.M_no_prev[, idx_stg3])
  pop_total_stg4 <- sum(m.M_no_prev[, idx_stg4])
  
  ## total de pacientes con cáncer clínico
  pop_total_cancer <- sum(m.M_no_prev[, idx_cancer])
  
  ## proporciones agregadas (1 valor por estadio)
  PropSick_estadio1_agg <- pop_total_stg1 / pop_total_cancer
  PropSick_estadio2_agg <- pop_total_stg2 / pop_total_cancer
  PropSick_estadio3_agg <- pop_total_stg3 / pop_total_cancer
  PropSick_estadio4_agg <- pop_total_stg4 / pop_total_cancer
  
  ## 7.1 target: distribucion estadios al diagnóstico
  v_ciclos <- 1:n.t
  
  v_diag_stg1 <- m.M_no_prev[v_ciclos, "1preclinical"] * a.P_noprev["1preclinical", "1clinicala", v_ciclos]
  v_diag_stg2 <- m.M_no_prev[v_ciclos, "2preclinical"] * a.P_noprev["2preclinical", "2clinicala", v_ciclos]
  v_diag_stg3 <- m.M_no_prev[v_ciclos, "3preclinical"] * a.P_noprev["3preclinical", "3clinicala", v_ciclos]
  v_diag_stg4 <- m.M_no_prev[v_ciclos, "4preclinical"] * a.P_noprev["4preclinical", "4clinicala", v_ciclos]

  
  ## la suma de los 4 estadios debe ser igual a los casos nuevos de cada ciclo
  stopifnot(isTRUE(all.equal(unname(v_diag_stg1 + v_diag_stg2 + v_diag_stg3 + v_diag_stg4),
                             v_new_cases[v_ciclos])))
  
  ## diagnósticos desde los 30 años.
  ciclos_diag <- 13:n.t
  
  diag_tot_stg1 <- sum(v_diag_stg1[ciclos_diag])
  diag_tot_stg2 <- sum(v_diag_stg2[ciclos_diag])
  diag_tot_stg3 <- sum(v_diag_stg3[ciclos_diag])
  diag_tot_stg4 <- sum(v_diag_stg4[ciclos_diag])
  
  ## total de diagnósticos en ciclos de interes
  diag_tot <- diag_tot_stg1 + diag_tot_stg2 + diag_tot_stg3 + diag_tot_stg4
  
  ## proporción de diagnósticos en cada estadio  (TARGET 3)
  PropDiag_estadio1 <- diag_tot_stg1 / diag_tot
  PropDiag_estadio2 <- diag_tot_stg2 / diag_tot
  PropDiag_estadio3 <- diag_tot_stg3 / diag_tot
  PropDiag_estadio4 <- diag_tot_stg4 / diag_tot
    
  ## 8. return results in a list
  out <- list(
    incidencia_edad     = v_tasa_inc_modelada,   ## TARGET: 6 tramos, pop viva  -> lst_targets$incidencia_edad
    incidencia_nacional = tasa_inc_nacional,     ## TARGET: ASR segi, pop viva  -> lst_targets$incidencia
    mortalidad          = v_tasa_mort_modelada,
    mortalidad_ciclo    = v_tasa_mort_ciclo,
    incidencia_ciclo    = v_tasa_inc_ciclo,
    
    ## tasas descriptivas sobre población en riesgo (No son target)
    incidencia_edad_riesgo     = v_tasa_inc_riesgo,
    incidencia_nacional_riesgo = tasa_inc_nacional_riesgo,
    incidencia_ciclo_riesgo    = v_tasa_inc_ciclo_riesgo,
    
    
    ## proporciones de estadios neoplásicos
    
    # Descriptivo: proporciones de Estadio ciclo a ciclo
    PropSick_estadio1 = v_prop_stg1,
    PropSick_estadio2 = v_prop_stg2,
    PropSick_estadio3 = v_prop_stg3,
    PropSick_estadio4 = v_prop_stg4,
    
    # Descriptivo: proporciones de Estadio Agregadas 
    PropSick_estadio1_agg = PropSick_estadio1_agg,
    PropSick_estadio2_agg = PropSick_estadio2_agg,
    PropSick_estadio3_agg = PropSick_estadio3_agg,
    PropSick_estadio4_agg = PropSick_estadio4_agg,
    
    ## TARGET 3 -> distribución de estadio al diagnóstico (30-99 años)
    PropDiag_estadio1 = PropDiag_estadio1,
    PropDiag_estadio2 = PropDiag_estadio2,
    PropDiag_estadio3 = PropDiag_estadio3,
    PropDiag_estadio4 = PropDiag_estadio4,
    
    ## diagnósticos por estadio en cada ciclo (82 valores, descriptivo)
    diag_estadio1 = unname(v_diag_stg1),
    diag_estadio2 = unname(v_diag_stg2),
    diag_estadio3 = unname(v_diag_stg3),
    diag_estadio4 = unname(v_diag_stg4),
    
    
    # Prevalencias Preneoplásicas 
    
    ## normal-gastritis (hp+ y hp-)
    PrevNG_combined = unname(sapply(l_filas_tramos_NG, function(filas) 
      prev_marginal(m.M_no_prev, 
                    c(idx_Nn, idx_Np, idx_Gn, idx_Gp),  ## numerador
                    idx_alive,                          ## denominador
                    filas))),                           ## filas del tramo
    
    PrevNG_combined_punt = calc_prev(m.M_no_prev, c(idx_Nn, idx_Np, idx_Gn, idx_Gp), ciclos_eval),
    
    ## metaplasia (hp+ y hp-)
    PrevI_combined  = prev_marginal(m.M_no_prev, 
                                    c(idx_In, idx_Ip), ## numerador
                                    idx_alive,         ## denominador
                                    rows_prev),        ## filas donde contar (años)
    
    ## displasia (hp+ y hp-)
    PrevD_combined  = prev_marginal(m.M_no_prev, c(idx_Dn, idx_Dp), idx_alive, rows_prev),
    
    ## atrofia hp+
    PrevA_hpp       = prev_marginal(m.M_no_prev, idx_Ap, idx_hpp, rows_prev),
    
    ## atrofia hp-
    PrevA_hpn       = prev_marginal(m.M_no_prev, idx_An, idx_hpn, rows_prev)
  
    )
  
  ## cohort trace, con argumento solo si se pide explícitamente (no usar al calibrar)
  if (return_trace) out$m.M_no_prev <- m.M_no_prev
  
  return(out)
  
  
  }  ## end as.list(), {}}
  
 )   ## end with() function
  
}    ## end runs_cg_markov() function