# 0. Libraries ------------------------------------------------------------
.libPaths(c("/projects/above_gedi/users/pburns/envs/biodiv_mod/lib/R/library","/scratch/pb463/R_libs"))

library(data.table)
library(randomForest)
library(prg) # from devtools
library(PresenceAbsence)
library(mRMRe)
library(Boruta)
library(collinear)
library(vip)
# reprtree funcs for converting RF model to text file
source('/home/pb463/scripts/repos/sunda_sdm_gedi/R/reprtree_update.R')
#library(ranger)
#library(caret)




# 1. Arguments ------------------------------------------------------------
run_local_test <- FALSE
if (run_local_test){
  # testing vals
  run_id <- 1
  region <- 'global'
  # region extent
  # currently not used since site coordinates have been removed
  #region_ext <- "108.44,119.54,-4.1,7.05" #xmin, xmax, ymin, ymax
  in_tab <- "/scratch/pb463/spec_pred_ext/borneo_pa_09072023_predext_stdz_sunda_clouded_leopard.csv"
  gedi_scen <- 'base'
  fs_method <- 'gvifcdr' 
  seed <- 123
  cores <- 4
  out_dir <- "/scratch/pb463/test/"
} else {
  
 # inputs
  args <- commandArgs(TRUE)
  # run identifier
  run_id <- as.numeric(args[1])
  
  # region
  region <- args[2]
  
  # region extent
  # currently not used since site coordinates have been removed
  region_ext <- args[3]
  
  # input table 
  in_tab <- args[4]
  
  # GEDI scenario
  gedi_scen <- args[5]
  
  # feature selection method
  fs_method <- args[6]

  # RF class balance
  rf_class_bal <- args[7]
  
  # random seed
  seed <- args[8]
  
  # cores
  cores <- args[9]
  
  # output base directory
  out_dir <- args[10]
}

# other hardcoded params

# the preference order metric when examining pairwise correlations, probably better to use eval metrics associated with training
# many options, for example: 'rf_aucroc', 'rfglmmn_aucroc_tr', 'rfglmmn_aucroc', 'rf_tss', 'glm_tss', 'glm_aucroc'
po_met <- 'aic_tss_sc_mn' 

# whether or not to use focal SD metrics
use_focsd <- TRUE

# whether or not to only use nominal scales (i.e. smallest) per predictor
nom_scales_only <- FALSE

# 2. Load Data ------------------------------------------------------------
# Load predictor extraction table
dt_o <- fread(in_tab)
f_split <- strsplit(x = tools::file_path_sans_ext(basename(in_tab)), split = "_")[[1]]
spec <- paste(f_split[6:length(f_split)], collapse="_")
pred_var_names <- colnames(dt_o)[grepl(pattern = "scl", x = colnames(dt_o), fixed = TRUE)]


out_base_name <- paste0("r", sprintf("%06d", run_id), "_", region, "_", gsub(pattern = "_", replacement = "-", x = spec), 
                        "_", gedi_scen, "_", fs_method, "_", rf_class_bal, "_", seed)

out_base_name_simp <- paste0(gsub(pattern = "_", replacement = "-", x = spec), 
                             "_", gedi_scen, "_", fs_method, "_", rf_class_bal, "_", seed)

# add x/y coordinates from .geo column
# dt_o <- dt_o[, c("delme", "xy") := tstrsplit(x = .geo, ":[", fixed=TRUE)]
# dt_o <- dt_o[, c("x","y"):=tstrsplit(xy, ",", fixed=TRUE)]
# dt_o$x <- as.numeric(dt_o$x)
# dt_o$y = as.numeric(gsub(pattern = "]|}", replacement = "", x = dt_o$y))
# dt_o <- dt_o[, c("delme", "xy") := NULL]

# remove unnecessary columns
dt_o <- dt_o[, c("system:index", ".geo") := NULL]



# 3. Preproc Data ------------------------------------------------------------
## get predictor variable names for the scenario 
# several options (some with two names per option)
# base = climate, disturbance, geomorph, human, productivity
# base_lsccdc = base + Landsat CCDC spectral indices
# base_gedicc = base + GEDI fusion veg. structure 
# ogedicc = only GEDI fusion veg. structure
# ogedikr = only GEDI kriged veg. structure
# hind_wgedicc (ghs) = high res. suitable for hindcasting (no climate, Landsat CCDC, or productivity)
# base_gedicc_msr = base + GEDI fusion veg. structure limited to focal means less than or equal to 1200 m
# full = everything except GEDI kr
if(gedi_scen == "nogedi" || gedi_scen == "base_lsccdc"){
  cat("Running base model, adding LS CCDC... \n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gkr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_ggr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gcc*"),pred_var_names)]
  
} else if (gedi_scen == "wgedicc" || gedi_scen == "base_gedicc"){
  cat("Running base model, adding GEDI CCDC fusion.. \n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gkr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_ggr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("comp*"),pred_var_names)]
  
} else if (gedi_scen == "wgedikr" || gedi_scen == "base_gedikr"){
  cat("Running base model, adding LS CCDC and GEDI kriged... \n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gcc*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_ggr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("comp*"),pred_var_names)]
  
} else if (gedi_scen == "wgedigr" || gedi_scen == "base_gedigr"){
  cat("Running base model, adding LS CCDC and GEDI gridded... \n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gcc*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gkr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("comp*"),pred_var_names)]
  
} else if (gedi_scen == "ogedicc"){
  cat("Running only GEDI CCDC fusion... \n")
  pred_var_names <- pred_var_names[grepl(glob2rx("*_gcc*"),pred_var_names)]
  
} else if (gedi_scen == "ogedikr"){
  cat("Running only GEDI kriged... \n")
  pred_var_names <- pred_var_names[grepl(glob2rx("*_gkr*"),pred_var_names)]
  
} else if (gedi_scen == "ogedigr"){
  cat("Running only GEDI gridded... \n")
  pred_var_names <- pred_var_names[grepl(glob2rx("*_ggr*"),pred_var_names)]
  
} else if (gedi_scen == "ghs" || gedi_scen == "hind_wgedicc"){
  search_str <- c('struct_anydeg', 'struct_anyloss', 'struct_recloss', 
                  'struct_gccrh95', 'struct_gccrh50', 'struct_gccpavd05', 'struct_gccpai', 'struct_gcccover', 'struct_gccfhd', 'struct_gccnummodes', 'struct_gccagbd', 
                  'hum_ghm', 'hum_pop', 'hum_gdp', 
                  'geo_swater', 'geo_elev', 'geo_slope', 'geo_lfpk', 'geo_lfsup', 'geo_lfslo', 'geo_lfval', 'geo_twi', 'geo_hand', 'geo_mtpi')
  p_list <- list()
  for (i in 1:length(search_str)){
    p_list[[i]] <- pred_var_names[grepl(glob2rx(pattern = paste0("*", search_str[i], "*")), pred_var_names)]
  }
  pred_var_names <- unlist(p_list)
  cat("Running with hindcast predictors (included GEDI CCDC fusion)... \n")
} else if (gedi_scen == 'base'){
  cat("Running base model... \n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gkr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_ggr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gcc*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("comp*"),pred_var_names)]
} else if (gedi_scen == 'base_gedicc_msr'){
  cat("Running base model, adding GEDI CCDC fusion only at fine to moderate spatial res... \n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gkr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_ggr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("comp*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*scl2400m"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*scl4800m"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*scl9600m"),pred_var_names)]
} else if (gedi_scen == 'full'){
  cat("Running full model (no GEDI kriging or gridding)... \n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_gkr*"),pred_var_names)]
  pred_var_names <- pred_var_names[!grepl(glob2rx("*_ggr*"),pred_var_names)]
} else {
  cat("gedi_scen argument is not valid...\n")
}

# whether or not to include focal SD predictors
if (use_focsd){
  cat("...including focal SD predictors...\n")
} else {
  cat("...excluding focal SD predictors...\n")
  pred_var_names <- pred_var_names[!grepl(glob2rx("*focsd*"),pred_var_names)]
}

# whether or not to only use nominal scales (i.e. smallest) per predictor
if (nom_scales_only){
  pred_var_names_nom <- c("clim_aridi_focmn_scl1200m", "clim_ccmean_focmn_scl1200m", "clim_hursmn_focmn_scl1200m", "clim_itherm_focmn_scl1200m", "clim_lstdiff_focmn_scl1200m", "clim_lstnimn_focmn_scl1200m", "clim_map_focmn_scl1200m", "clim_mat_focmn_scl1200m", "clim_petmn_focmn_scl1200m", "clim_pseas_focmn_scl1200m", "clim_rsdsmn_focmn_scl1200m", "clim_vpdmn_focmn_scl1200m", "clim_vpdrng_focmn_scl1200m", "clim_windmn_focmn_scl1200m", 
  "comp_nbr_nofoc_scl90m", "comp_ndmi_nofoc_scl90m", "comp_ndvi_nofoc_scl90m", "comp_svvi_nofoc_scl90m", 
  "ctns", "geo_clay_focmn_scl1200m", "geo_elev_nofoc_scl90m", "geo_hand_nofoc_scl90m", "geo_lfpk_focmn_scl150m", "geo_lfslo_focmn_scl150m", "geo_lfsup_focmn_scl150m", "geo_lfval_focmn_scl150m", "geo_mtpi_nofoc_scl90m", "geo_nitrogen_focmn_scl1200m", "geo_phh2o_focmn_scl1200m", "geo_sand_focmn_scl1200m", "geo_slope_nofoc_scl90m", "geo_soc_focmn_scl1200m", "geo_swater_focmn_scl150m", "geo_twi_nofoc_scl90m", 
  "hum_crops_focmn_scl150m", "hum_gdp_focmn_scl1200m", "hum_ghm_focmn_scl1200m", "hum_huntp_focmn_scl1200m", "hum_pa_focmn_scl1200m", "hum_palm_focmn_scl150m", 
  "prod_gppdhicum_focmn_scl1200m", "prod_gppdhimin_focmn_scl1200m", "prod_gppdhiseas_focmn_scl1200m", "prod_nppchels_focmn_scl1200m", 
  "struct_anydeg_focmn_scl150m", "struct_anyloss_focmn_scl150m", "struct_recloss_focmn_scl150m",
  "struct_gccagbd_nofoc_scl90m", "struct_gcccover_nofoc_scl90m", "struct_gccfhd_nofoc_scl90m", "struct_gccnummodes_nofoc_scl90m", "struct_gccpai_nofoc_scl90m", "struct_gccpavd05_nofoc_scl90m", "struct_gccrh50_nofoc_scl90m", "struct_gccrh95_nofoc_scl90m",
  "struct_gkragbd_focmn_scl1200m", "struct_gkrcover_focmn_scl1200m", "struct_gkrfhd_focmn_scl1200m", "struct_gkrnummodes_focmn_scl1200m", "struct_gkrpai_focmn_scl1200m", "struct_gkrpavd05_focmn_scl1200m", "struct_gkrrh50_focmn_scl1200m", "struct_gkrrh95_focmn_scl1200m",
  "struct_ggragbd_focmn_scl1200m", "struct_ggrcover_focmn_scl1200m", "struct_ggrfhd_focmn_scl1200m", "struct_ggrnmodes_focmn_scl1200m", "struct_ggrpai_focmn_scl1200m", "struct_ggrpavd05_focmn_scl1200m", "struct_ggrrh50_focmn_scl1200m", "struct_ggrrh95_focmn_scl1200m")
  pred_var_names <- pred_var_names[pred_var_names%in%pred_var_names_nom] 
}

# filter by region
# currently not used since site coordinates have been removed
# xmin <- as.numeric(strsplit(x = region_ext, split = ",", fixed = TRUE)[[1]][1])
# xmax <- as.numeric(strsplit(x = region_ext, split = ",", fixed = TRUE)[[1]][2])
# ymin <- as.numeric(strsplit(x = region_ext, split = ",", fixed = TRUE)[[1]][3])
# ymax <- as.numeric(strsplit(x = region_ext, split = ",", fixed = TRUE)[[1]][4])
# dt_s <- dt_o[x >= xmin & x <= xmax & y >= ymin & y <= ymax]

# omit NAs
dt_s <- na.omit(dt_o)

# make sure to include camera trap nights
pred_var_names <- c(pred_var_names, "ctns")
cat("There are", length(pred_var_names), "possible predictors \n")

# make sure column names are correct
# renaming ctns to match other naming convention
rseed<- as.character(paste0("gee_rseed",seed)) # name of the column for gee seed
keep_cols <- c(spec, pred_var_names, rseed, "location_id")
dt_s <- dt_s[, ..keep_cols]
colnames(dt_s) <- c("pa", pred_var_names, rseed, "location_id")

# check that there are enough presences for modeling
all_p_n <- nrow(dt_s[pa==1])
all_a_n <- nrow(dt_s[pa==0])

if (all_p_n < 20){
  cat("Not enough presences to run models. Quitting...")
  quit(save = "no", status = 0)
}


### Sampling for training and validation 25:75% -->
cat("Training / testing split \n")
# order by rseed
dt_s <- dt_s[order(get(rseed))]
#TRAINING dataset 75% pres and 75% abs BASED ON COLUMN SEED
#Presence
dt_t_p <- dt_s[pa==1][1:round(nrow(dt_s[pa==1])*0.75)]
#Absence
dt_t_a <- dt_s[pa==0][1:round(nrow(dt_s[pa==0])*0.75)]
#Combine
dt_t <- rbind(dt_t_p, dt_t_a)
rm(dt_t_p, dt_t_a)
cat(nrow(dt_t), "rows for training (75%) \n")

#VALIDATION dataset ~ 25% pres and 25% abs
dt_v <- dt_s[!(location_id %in% dt_t$location_id)] #delete rows used for training
dt_s <- dt_s[, c(rseed) := NULL]
dt_t <- dt_t[, c(rseed) := NULL]
dt_v <- dt_v[, c(rseed) := NULL]
cat(nrow(dt_v), "rows for testing (25%) \n")
cat("\n")

# check to see if there are any zero variance predictors in the training DT
var <- sapply(dt_t[,..pred_var_names], var)
zv_preds <- names(var[var<0.001])
cat("Removing near-zero variance predictors: ", zv_preds, "\n")
pred_var_names <- pred_var_names[!(pred_var_names %in% zv_preds)]
cat("\n")
# also check that each predictor has at least 3 unique values
luv <- sapply(dt_t[,..pred_var_names], function(x){length(unique(x))})
luv_preds <- names(luv[luv<=2])
luv_p <- sapply(dt_t[pa==1,..pred_var_names], function(x){length(unique(x))})
luv_p_preds <- names(luv_p[luv_p<=2])
cat("Removing predictors with less than or equal to 2 unique values: ", luv_preds, luv_p_preds, "\n")
pred_var_names <- pred_var_names[!((pred_var_names %in% luv_preds) | (pred_var_names %in% luv_p_preds))]
cat("\n")
keep_cols <- c("pa", pred_var_names, "location_id")
dt_t <- dt_t[, ..keep_cols]
dt_v <- dt_v[, ..keep_cols]

cat("There are", length(pred_var_names), "possible predictors after additional checks \n")


# function to fit and evaluate RF and GLM models
rf_glm_fit_eval <- function(y="pa", 
                            preds=NULL, 
                            dt_train=NULL, 
                            dt_val=NULL, 
                            #opt_met='ROC', 
                            mtry1 = FALSE,
                            ntrees = 100,
                            nodesize = 5,
                            rf_imp = FALSE,
                            save_rf_mod = FALSE,
                            #save_rfrang_mod = FALSE,
                            save_glm_mod = FALSE,
                            save_tr_fmt = FALSE){
  
  mod_cols <- c(y, preds)
  dt_m_t <- dt_train[ , ..mod_cols]
  colnames(dt_m_t) <- c("pa", preds)
  dt_m_v <- dt_val[ , ..mod_cols]
  colnames(dt_m_v) <- c("pa", preds)
  
  # get the observed validation values
  tr_obs <- as.numeric(dt_m_t$pa)
  te_obs <- as.numeric(dt_m_v$pa)
  
  # convert response column to factor in order to get classification probabilities
  #dt_m_t <- dt_m_t[,pa := ifelse(pa==1,"p","a")]
  dt_m_t$pa <- as.factor(x = dt_m_t$pa)
  #dt_m_v <- dt_m_v[,pa := ifelse(pa==1,"p","a")]
  dt_m_v$pa <- as.factor(x = dt_m_v$pa)
  
  
  # caret training specifications
  # eval_summ  <- function(data, lev = NULL, model = NULL){
  #   a1 <- defaultSummary(data, lev, model)
  #   b1 <- twoClassSummary(data, lev, model)
  #   prg_curve <- prg::create_prg_curve(labels = as.numeric(data$obs)-1, pos_scores = data$p)
  #   aucprg <- prg::calc_auprg(prg_curve)
  #   names(aucprg) <- 'auc_prg'
  #   tss <- b1['Sens'] + b1['Spec'] - 1
  #   names(tss) <- 'TSS'
  #   out <- c(a1, b1, aucprg, tss)
  #   out}
  # 
  # fitControl <- trainControl(
  #   method = "cv",
  #   number = 4,
  #   classProbs = TRUE,
  #   summaryFunction = eval_summ)
  
  # mtry specification
  if (mtry1){
    mtry <- 1
  } else {
    mtry_a <- max(c(1,round(length(preds)/3)))
    mtry_b <- max(c(1,round(sqrt(length(preds)))))
    mtrys <- sort(unique(c(mtry_a, mtry_b)))
    mtry <- mtry_b
  }
  
  if (rf_imp){
    rang_imp <- 'permutation'
  } else {
    rang_imp <- 'none'
  }
  # Ranger RF
  # tgrid_rang <- expand.grid(
  #   .mtry = mtrys,
  #   .splitrule = c("gini", "hellinger"),
  #   .min.node.size = c(2,4,8)
  # )
  
  # tr_fit_rang <- train(pa ~ ., data = dt_m_t, 
  #                      method = "ranger", 
  #                      trControl = fitControl,
  #                      tuneGrid = tgrid_rang,
  #                      metric = opt_met,
  #                      num.trees = ntrees,
  #                      importance = rang_imp,
  #                      probability = TRUE,
  #                      verbose = FALSE)
  
  # Original RF
  # tgrid_rf <- expand.grid(.mtry = mtrys)
  # tr_fit_rf <- train(pa ~ ., data = dt_m_t,
  #                    method = 'rf',
  #                    trControl = fitControl,
  #                    tuneGrid = tgrid_rf,
  #                    ntree=ntrees,
  #                    importance = TRUE,
  #                    nodesize=2,
  #                    metric = opt_met,
  #                    verbose = FALSE)
  
  if (rf_class_bal == 'bal'){
  # class balance
  samp_n <- min(nrow(dt_m_t[pa=='1']), nrow(dt_m_t[pa=='0']))
  tr_fit_rfds <- randomForest(pa~., data = dt_m_t, importance = rf_imp,
                              ntree = ntrees, mtry = mtry, nodesize = nodesize,
                              sampsize = c("0" = samp_n, "1" = samp_n))
  } else {
  tr_fit_rfds <- randomForest(pa~., data = dt_m_t, importance = rf_imp,
                              ntree = ntrees, mtry = mtry, nodesize = nodesize)    
  }
  
  # GLM
  tr_fit_glm <- glm(pa ~ ., data=dt_m_t, family="binomial")
  
  # Make predictions on training and test sets
  #te_pred_rang <- predict(tr_fit_rang, newdata = dt_m_v, type = 'prob')[,'p']
  tr_pred_rf <- predict(tr_fit_rfds, newdata = dt_m_t, type = 'prob')[,'1']
  tr_pred_glm <- as.vector(predict(tr_fit_glm, newdata = dt_m_t, type = 'response'))
  te_pred_rf <- predict(tr_fit_rfds, newdata = dt_m_v, type = 'prob')[,'1']
  te_pred_glm <- as.vector(predict(tr_fit_glm, newdata = dt_m_v, type = 'response'))
  
  # calc AUPRG curve
  #prg_curve_rang <- create_prg_curve(labels = te_obs, pos_scores = te_pred_rang)
  #aucprg_rang <- calc_auprg(prg_curve_rang)
  tr_prg_curve_rf <- create_prg_curve(labels = tr_obs, pos_scores = tr_pred_rf)
  tr_aucprg_rf <- calc_auprg(tr_prg_curve_rf)
  te_prg_curve_rf <- create_prg_curve(labels = te_obs, pos_scores = te_pred_rf)
  te_aucprg_rf <- calc_auprg(te_prg_curve_rf)
  tr_prg_curve_glm <- create_prg_curve(labels = tr_obs, pos_scores = tr_pred_glm)
  tr_aucprg_glm <- calc_auprg(tr_prg_curve_glm)
  te_prg_curve_glm <- create_prg_curve(labels = te_obs, pos_scores = te_pred_glm)
  te_aucprg_glm <- calc_auprg(te_prg_curve_glm)
  
  # calc AUC and other metrics
  tr_id <- seq_len(length(tr_obs))
  te_id <- seq_len(length(te_obs))
  # ranger first
  #te_df_rang <- data.frame(id = id, obs = te_obs, pred = te_pred_rang)
  #acc_rang <- presence.absence.accuracy(te_df_rang, threshold=500, st.dev=FALSE)
  #maxacc_rang <- acc_rang[acc_rang[,'Kappa'] == max(acc_rang[,'Kappa']),] 
  #maxacc_rang <- maxacc_rang[1,2:7]
  #maxacc_rang$TSS <- maxacc_rang$sensitivity+maxacc_rang$specificity-1
  #names(maxacc_rang) <- c("rfrang_aucrocth", "rfrang_pcc", "rfrang_sens", "rfrang_spec", 
  #                        "rfrang_kappa", "rfrang_aucroc", "rfrang_tss")
  #maxacc_rang$rfrang_aucprg <- aucprg_rang
  #maxacc_rang$rfrang_rocprgmn <- mean(c(maxacc_rang$rfrang_aucroc, aucprg_rang))
  
  # RF, training and testing
  tr_df_rf <- data.frame(id = tr_id, obs = tr_obs, pred = tr_pred_rf)
  te_df_rf <- data.frame(id = te_id, obs = te_obs, pred = te_pred_rf)
  
  tr_acc_rf <- presence.absence.accuracy(tr_df_rf, threshold=200, st.dev=FALSE)
  tr_maxacc_rf <- tr_acc_rf[tr_acc_rf[,'Kappa'] == max(tr_acc_rf[,'Kappa']),] 
  tr_maxacc_rf <- tr_maxacc_rf[1,2:7]
  tr_maxacc_rf$TSS <- tr_maxacc_rf$sensitivity+tr_maxacc_rf$specificity-1
  names(tr_maxacc_rf) <- c("rf_aucrocth_tr", "rf_pcc_tr", "rf_sens_tr", "rf_spec_tr", 
                           "rf_kappa_tr", "rf_aucroc_tr", "rf_tss_tr")
  tr_maxacc_rf$rf_aucprg_tr <- tr_aucprg_rf
  tr_maxacc_rf$rf_rocprgmn_tr <- mean(c(tr_maxacc_rf$rf_aucroc_tr, tr_aucprg_rf))
  tr_maxacc_rf$rf_erate0_tr <- tr_fit_rfds$confusion[5]
  tr_maxacc_rf$rf_erate1_tr <- tr_fit_rfds$confusion[6]
  tr_maxacc_rf$rf_sens_troob <- tr_fit_rfds$confusion[4]/(tr_fit_rfds$confusion[4]+tr_fit_rfds$confusion[2])
  tr_maxacc_rf$rf_spec_troob <- tr_fit_rfds$confusion[1]/(tr_fit_rfds$confusion[1]+tr_fit_rfds$confusion[3])
  tr_maxacc_rf$rf_prec_troob <- tr_fit_rfds$confusion[4]/(tr_fit_rfds$confusion[4]+tr_fit_rfds$confusion[3])
  tr_maxacc_rf$rf_tss_troob <- (tr_maxacc_rf$rf_sens_troob + tr_maxacc_rf$rf_spec_troob - 1)
  
  te_acc_rf <- presence.absence.accuracy(te_df_rf, threshold=200, st.dev=FALSE)
  te_maxacc_rf <- te_acc_rf[te_acc_rf[,'Kappa'] == max(te_acc_rf[,'Kappa']),] 
  te_maxacc_rf <- te_maxacc_rf[1,2:7]
  te_maxacc_rf$TSS <- te_maxacc_rf$sensitivity+te_maxacc_rf$specificity-1
  names(te_maxacc_rf) <- c("rf_aucrocth", "rf_pcc", "rf_sens", "rf_spec", 
                           "rf_kappa", "rf_aucroc", "rf_tss")
  te_maxacc_rf$rf_aucprg <- te_aucprg_rf
  te_maxacc_rf$rf_rocprgmn <- mean(c(te_maxacc_rf$rf_aucroc, te_aucprg_rf))
  
  
  # GLM, training and testing
  tr_df_glm <- data.frame(id = tr_id, obs = tr_obs, pred = tr_pred_glm)
  tr_acc_glm <- presence.absence.accuracy(tr_df_glm, threshold=200, st.dev=FALSE)
  tr_maxacc_glm <- tr_acc_glm[tr_acc_glm[,'Kappa'] == max(tr_acc_glm[,'Kappa']),] 
  tr_maxacc_glm <- tr_maxacc_glm[1,2:7]
  tr_maxacc_glm$TSS <- tr_maxacc_glm$sensitivity+tr_maxacc_glm$specificity-1
  names(tr_maxacc_glm) <- c("glm_aucrocth_tr", "glm_pcc_tr", "glm_sens_tr", "glm_spec_tr", 
                            "glm_kappa_tr", "glm_aucroc_tr", "glm_tss_tr")
  tr_maxacc_glm$glm_aucprg_tr <- tr_aucprg_glm
  tr_maxacc_glm$glm_rocprgmn_tr <- mean(c(tr_maxacc_glm$glm_aucroc_tr, tr_aucprg_glm))
  tr_maxacc_glm$glm_aic_tr <- tr_fit_glm$aic
  
  te_df_glm <- data.frame(id = te_id, obs = te_obs, pred = te_pred_glm)
  te_acc_glm <- presence.absence.accuracy(te_df_glm, threshold=200, st.dev=FALSE)
  te_maxacc_glm <- te_acc_glm[te_acc_glm[,'Kappa'] == max(te_acc_glm[,'Kappa']),] 
  te_maxacc_glm <- te_maxacc_glm[1,2:7]
  te_maxacc_glm$TSS <- te_maxacc_glm$sensitivity+te_maxacc_glm$specificity-1
  names(te_maxacc_glm) <- c("glm_aucrocth", "glm_pcc", "glm_sens", "glm_spec", 
                            "glm_kappa", "glm_aucroc", "glm_tss")
  te_maxacc_glm$glm_aucprg <- te_aucprg_glm
  te_maxacc_glm$glm_rocprgmn <- mean(c(te_maxacc_glm$glm_aucroc, te_aucprg_glm))
  te_maxacc_glm$glm_aic <- tr_fit_glm$aic
  
  dt_eval <- as.data.table(cbind(tr_maxacc_rf, te_maxacc_rf, tr_maxacc_glm, te_maxacc_glm))
  # compute average RF + GLM aucroc and TSS
  dt_eval <- dt_eval[,rfglmmn_aucroc_tr := fifelse(!is.na(rf_aucroc_tr) & !is.na(glm_aucroc_tr), 
                                                   ((rf_aucroc_tr + glm_aucroc_tr) / 2),
                                                   rf_aucroc_tr)]
  dt_eval <- dt_eval[,rfglmmn_aucroc := fifelse(!is.na(rf_aucroc) & !is.na(glm_aucroc), 
                                                ((rf_aucroc + glm_aucroc) / 2),
                                                rf_aucroc)]
  
  dt_eval <- dt_eval[,rfglmmn_tss_tr := fifelse(!is.na(rf_tss_tr) & !is.na(glm_tss_tr), 
                                                ((rf_tss_tr + glm_tss_tr) / 2),
                                                rf_tss_tr)]
  dt_eval <- dt_eval[,rfglmmn_tss := fifelse(!is.na(rf_tss) & !is.na(glm_tss), 
                                             ((rf_tss + glm_tss) / 2),
                                             rf_tss)]
  
  dt_eval$preds <- paste0(preds, collapse=",")
  
  # if (save_rfrang_mod){
  #   rfrang_mod <- tr_fit_rang$finalModel
  # } else {
  #   rfrang_mod <- NULL
  # }
  
  if (save_rf_mod){
    rf_mod <- tr_fit_rfds
  } else {
    rf_mod <- NULL
  }
  
  if (save_glm_mod){
    glm_mod <- tr_fit_glm
  } else {
    glm_mod <- NULL
  }
  
  if (save_tr_fmt){
    dt_tr_fmt <- dt_m_t
  } else {
    dt_tr_fmt <- NULL
  }
  
  return(list(dt_eval, rf_mod, glm_mod, dt_tr_fmt))
}


# functions for predicting RF and GLM
pfun_rf_prob <- function(object, newdata) {  
  pred_p_rf <- predict(object, newdata = newdata, type = "prob")[, "1"]
  return(pred_p_rf)
}

pfun_glm_prob <- function(object, newdata) {  
  pred_p_glm<- as.vector(predict(object, newdata = newdata, type = 'response'))
  return(pred_p_glm)
}

# function to scale values between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}



# 4. Feature selection ------------------------------------------------------------
# Functions for selecting predictors 
# MRMR feature selection function
# returns the top n feature names
mrmr_fs <- function(dt_in = NULL,
                    preds_vec = NULL,
                    max_n_sel = 15){
  mod_cols <- c('pa', preds_vec)
  dd <- dt_in[,..mod_cols]
  dd$pa <- as.numeric(dd$pa)
  dm <- mRMR.data(data = dd)
  mrmr_e <- mRMR.ensemble(data = dm, target_indices = c(1), solution_count = 20, feature_count = max_n_sel)
  mrmr_e_all_vars <- as.vector(matrix(mrmr_e@filters$`1`))
  top_n_vars <- sort(table(mrmr_e_all_vars), decreasing=TRUE)[1:max_n_sel]
  preds_fs <- sort(mrmr_e@feature_names[as.numeric(names(top_n_vars))])
  
  cat("Filtered predictors using MRMR. Reduced from", length(preds_vec),
      "to", length(preds_fs), "predictors \n")
  cat(preds_fs,"\n")
  return(as.vector(preds_fs))
}

# use VIF with preference order of all variables at once (no groups)
avif_fs <- function(dt_in = dt_t,
                    preds_vec = NULL,
                    dt_uni = dt_u,
                    po_col = po_met,
                    max_n_sel = 10){
  
  # reduce correlation in predictors with cor_select()
  # avoid rare cases where predictors are perfectly correlated
  preds_cs <- cor_select(df = dt_in,
                         predictors = preds_vec,
                         preference_order = dt_uni[order(get(po_col), decreasing = TRUE), preds],
                         max_cor = 0.95)
  
  preds_vif_sel <- collinear::vif_select(df = dt_in[,..preds_cs], 
                                         predictors = preds_cs,
                                         preference_order = dt_uni[order(get(po_col), decreasing = TRUE), preds], 
                                         max_vif = 5)
  if (length(preds_vif_sel) <= max_n_sel){
    preds_fs <- preds_vif_sel
  } else {
    # reduce number of selected variables further using univariate RF importance and eval
    # get RF importance and then rank
    fit_mods <- rf_glm_fit_eval(y = "pa", preds = preds_vif_sel, dt_train = dt_in, dt_val = dt_in, 
                                rf_imp = FALSE,
                                save_rf_mod = TRUE, save_glm_mod = TRUE)
    
    dt_in$pa <- as.factor(dt_in$pa)
    # RF AUC perm imp
    dt_rfimp <- as.data.table(vi(fit_mods[[2]],
                                 method = "permute",
                                 train = dt_in,
                                 target = "pa",
                                 feature_names = preds_vif_sel,
                                 metric = "roc_auc",
                                 event_level = "second",
                                 pred_wrapper = pfun_rf_prob,
                                 nsim = 50,
                                 type = 'difference'))
    colnames(dt_rfimp) <- c("preds", "rfimp_mndauc", "rfimp_sddauc")
    dt_rfimp <- dt_rfimp[,rfimp_mndauc_adj:=fifelse(rfimp_mndauc<0,0,rfimp_mndauc)]
    
    # GLM AUC perm imp
    dt_glmimp <- as.data.table(vi(fit_mods[[3]],
                                  method = "permute",
                                  train = dt_in,
                                  target = "pa",
                                  feature_names = preds_vif_sel,
                                  metric = "roc_auc",
                                  event_level = "second",
                                  pred_wrapper = pfun_glm_prob,
                                  nsim = 50,
                                  type = 'difference'))
    colnames(dt_glmimp) <- c("preds", "glmimp_mndauc", "glmimp_sddauc")
    dt_glmimp <- dt_glmimp[,glmimp_mndauc_adj:=fifelse(glmimp_mndauc<0,0,glmimp_mndauc)]
    
    # join AUC perm imp tables
    dt_imp <- dt_glmimp[dt_rfimp, on = "preds"]
    
    # average the decrease in AUC from the two models
    # also scale values between 0 and 1
    dt_imp <- dt_imp[,rfglm_mndauc_adj := ((rfimp_mndauc_adj + glmimp_mndauc_adj)/2)]
    dt_imp <- dt_imp[,rfglm_mndauc_adj_sc := range01(rfglm_mndauc_adj)]
    dt_imp <- dt_imp[order(rfglm_mndauc_adj_sc, decreasing = TRUE)]
    dt_imp <- dt_imp[,i_rank:=.I]
    
    # get univariate eval metrics and then rank
    dt_u_p <- dt_uni[preds %in% preds_vif_sel]
    dt_u_p <- dt_u_p[,rfglmmn_aucroc_adj:=fifelse(rfglmmn_aucroc<0.5,0,rfglmmn_aucroc-0.5)]
    dt_u_p <- dt_u_p[,rfglmmn_aucroc_adj_sc := range01(rfglmmn_aucroc_adj)]
    dt_u_p <- dt_u_p[order(get(po_col), decreasing = TRUE)]
    dt_u_p <- dt_u_p[,e_rank:=.I]
    
    # join importance and eval, then average scaled importance and eval metrics
    dt_j <- dt_imp[dt_u_p, on = "preds"]
    dt_j <- dt_j[,ie_rank:=((i_rank + e_rank)/2)]
    dt_j <- dt_j[,ie_sc:=((rfglm_mndauc_adj_sc + get(po_col))/2)]
    
    # order by average RF importance + eval scaled values, then by univariate eval
    dt_j <- dt_j[order(ie_sc, decreasing = TRUE)]
    preds_fs <- dt_j[1:max_n_sel, preds]
  }
  return(preds_fs)
}

# use VIF by variable group with preference order
# use VIF with preference order of all variables at once (no groups)
gvif_fs <- function(dt_in = dt_t,
                    preds_vec = NULL,
                    dt_uni = dt_u,
                    po_col = po_met,
                    max_n_sel = 10){
  
  # reduce correlation in predictors with cor_select()
  # avoid rare cases where predictors are perfectly correlated
  preds_cs <- cor_select(df = dt_in,
                         predictors = preds_vec,
                         preference_order = dt_uni[order(get(po_col), decreasing = TRUE), preds],
                         max_cor = 0.95)
  
  # split up variable names
  preds_sp <- tstrsplit(preds_cs[preds_cs != "ctns"], "_", fixed = TRUE)
  groups <- sort(unique(preds_sp[[1]]))
  vars <- sort(unique(preds_sp[[2]]))
  
  # first find the best scale+focal stat per variable, be very restrictive
  var_sel_g_comb <- c()
  for (g in groups){
    var_sel_v_comb <- c()
    g_vars <- unique(tstrsplit(x = preds_cs[grepl(pattern = paste0(g,"_"), x = preds_cs)], "_", fixed = TRUE)[[2]])
    for (v in g_vars){
      vs_in_g <- preds_cs[grepl(pattern = v, x = preds_cs)]
      if (length(vs_in_g) > 1){
        var_sel_v <- collinear::vif_select(df = dt_in[,..vs_in_g], 
                                           predictors = vs_in_g,
                                           preference_order = dt_uni[order(get(po_col), decreasing = TRUE), preds], 
                                           max_vif = 2.5)
        var_sel_v_comb <- c(var_sel_v_comb, var_sel_v)
      } else if (length(vs_in_g) == 1){
        # need this scenario in case there's only one scale/focstat per variable
        var_sel_v_comb <- c(var_sel_v_comb, vs_in_g)
      }
    }
    
    # now find the best variables per group, be moderately restrictive
    vs_by_g <- sort(var_sel_v_comb)
    if (length(vs_by_g) > 1){
      var_sel_g <- collinear::vif_select(df = dt_in[,..vs_by_g], 
                                         predictors = vs_by_g,
                                         preference_order = dt_uni[order(get(po_col), decreasing = TRUE), preds], 
                                         max_vif = 5)
      var_sel_g_comb <- c(var_sel_g_comb, var_sel_g)
    } else if (length(vs_by_g) == 1){
      var_sel_g_comb <- c(var_sel_g_comb, vs_by_g)
    }
  }
  
  # now find the best variables of the combined group results
  # need to add ctns back in
  gvif_keep <- sort(c('ctns',var_sel_g_comb))
  preds_vif_sel <- collinear::vif_select(df = dt_in[,..gvif_keep], 
                                         predictors = gvif_keep,
                                         preference_order = dt_uni[order(get(po_col), decreasing = TRUE), preds], 
                                         max_vif = 10)
  
  if (length(preds_vif_sel) <= max_n_sel){
    preds_fs <- preds_vif_sel
  } else {
    # reduce number of selected variables further using univariate RF importance and eval
    # get RF importance and then rank
    fit_mods <- rf_glm_fit_eval(y = "pa", preds = preds_vif_sel, dt_train = dt_in, dt_val = dt_in,
                                rf_imp = FALSE,
                                save_rf_mod = TRUE, save_glm_mod = TRUE)
    
    dt_in$pa <- as.factor(dt_in$pa)
    # RF AUC perm imp
    dt_rfimp <- as.data.table(vi(fit_mods[[2]],
                                 method = "permute",
                                 train = dt_in,
                                 target = "pa",
                                 feature_names = preds_vif_sel,
                                 metric = "roc_auc",
                                 event_level = "second",
                                 pred_wrapper = pfun_rf_prob,
                                 nsim = 50,
                                 type = 'difference'))
    colnames(dt_rfimp) <- c("preds", "rfimp_mndauc", "rfimp_sddauc")
    dt_rfimp <- dt_rfimp[,rfimp_mndauc_adj:=fifelse(rfimp_mndauc<0,0,rfimp_mndauc)]
    
    # GLM AUC perm imp
    dt_glmimp <- as.data.table(vi(fit_mods[[3]],
                                  method = "permute",
                                  train = dt_in,
                                  target = "pa",
                                  feature_names = preds_vif_sel,
                                  metric = "roc_auc",
                                  event_level = "second",
                                  pred_wrapper = pfun_glm_prob,
                                  nsim = 50,
                                  type = 'difference'))
    colnames(dt_glmimp) <- c("preds", "glmimp_mndauc", "glmimp_sddauc")
    dt_glmimp <- dt_glmimp[,glmimp_mndauc_adj:=fifelse(glmimp_mndauc<0,0,glmimp_mndauc)]
    
    # join AUC perm imp tables
    dt_imp <- dt_glmimp[dt_rfimp, on = "preds"]
    
    # average the decrease in AUC from the two models
    # also scale values between 0 and 1
    dt_imp <- dt_imp[,rfglm_mndauc_adj := ((rfimp_mndauc_adj + glmimp_mndauc_adj)/2)]
    dt_imp <- dt_imp[,rfglm_mndauc_adj_sc := range01(rfglm_mndauc_adj)]
    dt_imp <- dt_imp[order(rfglm_mndauc_adj_sc, decreasing = TRUE)]
    dt_imp <- dt_imp[,i_rank:=.I]
    
    # get univariate eval metrics and then rank
    dt_u_p <- dt_uni[preds %in% preds_vif_sel]
    dt_u_p <- dt_u_p[,rfglmmn_aucroc_adj:=fifelse(rfglmmn_aucroc<0.5,0,rfglmmn_aucroc-0.5)]
    dt_u_p <- dt_u_p[,rfglmmn_aucroc_adj_sc := range01(rfglmmn_aucroc_adj)]
    dt_u_p <- dt_u_p[order(get(po_col), decreasing = TRUE)]
    dt_u_p <- dt_u_p[,e_rank:=.I]
    
    # join importance and eval, then average scaled importance and eval metrics
    dt_j <- dt_imp[dt_u_p, on = "preds"]
    dt_j <- dt_j[,ie_rank:=((i_rank + e_rank)/2)]
    dt_j <- dt_j[,ie_sc:=((rfglm_mndauc_adj_sc + get(po_col))/2)]
    
    # order by average RF importance + eval scaled values, then by univariate eval
    dt_j <- dt_j[order(ie_sc, decreasing = TRUE)]
    preds_fs <- dt_j[1:max_n_sel, preds]
  }
  return(preds_fs)
}

# RF+GLM custom dredge
cdredge_fs <- function(dt_train = dt_t,
                       dt_val = dt_v,
                       y = "pa",
                       preds_vec = NULL,
                       max_n_preds = 10,
                       pred_sel_m = "freq" # options are "best" or "freq"
){
  mod_cols <- c(y, preds_vec)
  dt_m_t <- dt_train[ , ..mod_cols]
  colnames(dt_m_t) <- c("pa", preds_vec)
  dt_m_v <- dt_val[ , ..mod_cols]
  colnames(dt_m_v) <- c("pa", preds_vec)
  
  # get the observed validation values
  te_obs <- as.numeric(dt_m_v$pa)
  id <- seq_len(length(te_obs))
  
  # convert response column to factor in order to get classification probabilities
  #dt_m_t <- dt_m_t[,pa := ifelse(pa==1,"p","a")]
  dt_m_t$pa <- as.factor(x = dt_m_t$pa)
  #dt_m_v <- dt_m_v[,pa := ifelse(pa==1,"p","a")]
  dt_m_v$pa <- as.factor(x = dt_m_v$pa)
  
  # make an empty list to store dredge results
  dredge_list <- list()
  i <- 1
  system.time(for (v in 2:max_n_preds){
    cat("Dredging all combinations of models with", v, "predictors... \n")
    combos <- combn(x = seq(1:length(preds_vec)), m = v)
    cat(ncol(combos), "combinations \n")
    for (n in 1:ncol(combos)){
      preds_idx <- combos[,n]
      preds_idx_names <- preds_vec[preds_idx]
      mod_cols <- c("pa", preds_idx_names)
      samp_n <- min(nrow(dt_m_t[pa=='1']), nrow(dt_m_t[pa=='0']))
      if (rf_class_bal == 'bal'){
      rf_fit <- randomForest::randomForest(pa~., data = dt_m_t[,..mod_cols],
                                           ntree=100, mtry = round(sqrt(length(preds_idx_names))),
                                           nodesize=5,
                                           sampsize = c("0" = samp_n, "1" = samp_n),
                                           importance = FALSE,
                                           proximity = FALSE)
      } else {
      rf_fit <- randomForest::randomForest(pa~., data = dt_m_t[,..mod_cols],
                                           ntree=100, mtry = round(sqrt(length(preds_idx_names))),
                                           nodesize=5,
                                           importance = FALSE,
                                           proximity = FALSE)
      }
      glm_fit <- glm(pa~., data = dt_m_t[,..mod_cols], family = "binomial")
      
      #te_pred_rf <- predict(rf_fit, newdata = dt_m_v[,..mod_cols], type = 'prob')[,'1']
      te_pred_rf <- pfun_rf_prob(rf_fit, dt_m_v)
      #te_pred_glm <- as.vector(predict(glm_fit, newdata = dt_m_v[,..mod_cols], type = 'response'))
      te_pred_glm <- pfun_glm_prob(glm_fit, dt_m_v)
      
      dredge_list[[i]] <- data.table(i = i,
                                     n = n, 
                                     preds_idx = paste0(preds_idx, collapse = ","),
                                     rf_aucroc = PresenceAbsence::auc(data.table(id = id, obs = te_obs, pred = te_pred_rf), st.dev=FALSE),
                                     glm_aucroc = PresenceAbsence::auc(data.table(id = id, obs = te_obs, pred = te_pred_glm), st.dev=FALSE))
      i <- i+1
    }
    cat("\n")
  })
  
  dt_c <- data.table::rbindlist(dredge_list)
  
  # compute mean RF+GLM AUC and sort 
  dt_c <- dt_c[,rfglmmn_aucroc := ((rf_aucroc + glm_aucroc)/2)]
  dt_c <- dt_c[order(rfglmmn_aucroc, decreasing = TRUE)]
  
  # two option for getting the best variables index
  if (pred_sel_m == "best"){
    # 1. get the predictors associated with the highest AUC model
    vi <- as.integer(unlist(strsplit(x = dt_c[1,preds_idx], split = ",")))
  } else if (pred_sel_m == "freq"){
    # 2. get the most frequently used predictors from the top models
    auc_q99 <- as.vector(quantile(x = dt_c$rfglmmn_aucroc, probs = c(0.99)))
    n_mods <- nrow(dt_c[rfglmmn_aucroc >= auc_q99])
    keep_thresh <- 0.1
    p_list <- list()
    for (i in 1:n_mods){
      p_list[[i]] <- as.numeric(as.vector(strsplit(x = as.character(dt_c[i,3]), split = ",", fixed = TRUE)[[1]]))
    }
    freq_sort <- sort(table(unlist(p_list)), decreasing = TRUE)/n_mods
    freq_idx <- as.numeric(names(freq_sort[freq_sort >= keep_thresh]))
    freq_idx_len <- length(freq_idx)
    
    if (freq_idx_len > max_n_preds){
      vi <- freq_idx[1:max_n_preds]
    } else {
      vi <- freq_idx
    }
  } else {
    cat("Invalid pred_sel_m argument \n")
  }
  # return the predictors associated with the index
  return(preds_vec[vi])
}


# run univariate models if necessary
if (fs_method %in% c('avif', 'gvif', 'avifcdr', 'gvifcdr')){
  # Run univariate models
  uni_mod_res_list <- list()
  for (i in 1:length(pred_var_names)){
    pred_var_i <- pred_var_names[i]
    cat("Working on predictor", i, "-", pred_var_i, "\n")
    uni_mod_out <- rf_glm_fit_eval(y = "pa", preds = pred_var_i, 
                                   dt_train = dt_t, dt_val = dt_v, 
                                   mtry1 = TRUE, rf_imp = FALSE, ntrees = 100, 
                                   save_rf_mod = TRUE, save_glm_mod = FALSE, save_tr_fmt = FALSE)
    uni_mod_res_list[[i]] <- uni_mod_out[[1]]
  }
  # combine the univariate eval metrics from all preds into one data.table
  dt_u <- rbindlist(uni_mod_res_list)
  
  # compute scaled and inverted GLM AIC and RF error rate of class 1 for selecting the best preds
  dt_u <- dt_u[,glm_aic_tr_sc:=(1-range01(glm_aic_tr))]
  dt_u <- dt_u[,rf_erate1_tr_sc:=(1-range01(rf_erate1_tr))]
  dt_u <- dt_u[,aic_erate1_sc_mn:=(glm_aic_tr_sc + rf_erate1_tr_sc)/2]
  dt_u <- dt_u[,rf_tss_troob_sc:=range01(rf_tss_troob)]
  dt_u <- dt_u[,aic_tss_sc_mn:=(glm_aic_tr_sc + rf_tss_troob_sc)/2]
  
  # save the univariate eval metrics
  fwrite(x = dt_u, file = paste0(out_dir, 'uni/', out_base_name, '.csv'))
}


# run the specified feature selection
if (fs_method == 'mrmr'){
  # Only use MRMR for feature selection
  pred_var_names_fs = mrmr_fs(dt_in = dt_t, 
                              preds_vec = pred_var_names, 
                              max_n_sel = min(length(pred_var_names), 20))
  
} else if (fs_method == 'avif'){
  # Use all at once VIF 
  pred_var_names_fs = avif_fs(dt_in = dt_t, 
                              preds_vec = pred_var_names, 
                              dt_uni = dt_u, po_col = po_met, 
                              max_n_sel = 20)
  
} else if (fs_method == 'gvif'){
  # Use VIF by group 
  pred_var_names_fs = gvif_fs(dt_in = dt_t, 
                              preds_vec = pred_var_names, 
                              dt_uni = dt_u, po_col = po_met, 
                              max_n_sel = 20)
  
} else if (fs_method == 'mrmrcdr'){
  # Use MRMR with custom RF+GLM dredge
  pred_var_names_f = mrmr_fs(dt_in = dt_t, 
                             preds_vec = pred_var_names, 
                             max_n_sel = min(length(pred_var_names), 15))
  
  pred_var_names_fs <- cdredge_fs(dt_train = dt_t, dt_val = dt_v, y = "pa", 
                                  preds_vec = pred_var_names_f, 
                                  max_n_preds = min(c(10, length(pred_var_names_f))),
                                  pred_sel_m = "freq")
  
} else if (fs_method == 'avifcdr'){
  # Use all at once VIF with custom RF+GLM dredge
  pred_var_names_f = avif_fs(dt_in = dt_t, 
                             preds_vec = pred_var_names, dt_uni = dt_u, 
                             po_col = po_met, 
                             max_n_sel = 15)
  
  pred_var_names_fs <- cdredge_fs(dt_train = dt_t, dt_val = dt_v, y = "pa", 
                                  preds_vec = pred_var_names_f, 
                                  max_n_preds = min(c(10, length(pred_var_names_f))),
                                  pred_sel_m = "freq")
  
} else if (fs_method == 'gvifcdr'){
  # Use VIF by group with custom RF+GLM dredge
  pred_var_names_f = gvif_fs(dt_in = dt_t, 
                             preds_vec = pred_var_names, dt_uni = dt_u, 
                             po_col = po_met, 
                             max_n_sel = 15)
  
  pred_var_names_fs <- cdredge_fs(dt_train = dt_t, dt_val = dt_v, y = "pa", 
                                  preds_vec = pred_var_names_f, 
                                  max_n_preds = min(c(10, length(pred_var_names_f))),
                                  pred_sel_m = "freq")
}

# save the selected features 
fwrite(x = data.table(sel_pred_vars = pred_var_names_fs), 
       file = paste0(out_dir, 'fs/', out_base_name, '.csv'))



# 5. Fit, evaluate, and interpret the final models ------------------------------------------------------------
final_mods <- rf_glm_fit_eval(y = "pa", preds = pred_var_names_fs, 
                              dt_train = dt_t, dt_val = dt_v,
                              mtry1 = FALSE, ntrees = 200, 
                              rf_imp = TRUE, save_rf_mod = TRUE, save_glm_mod = TRUE, save_tr_fmt = TRUE)

# save model training and validation predicted probabilities
dt_t$tv <- "t"
dt_v$tv <- "v"
dt_tv <- rbind(dt_t, dt_v)
dt_p <- dt_tv[,..pred_var_names_fs]
dt_p$pa <- as.factor(dt_tv$pa) 
dt_p$tv <- dt_tv$tv
dt_p$location_id <- dt_tv$location_id
dt_p$rf_prob1 <- pfun_rf_prob(final_mods[[2]], newdata = dt_tv)
dt_p$glm_prob1 <- pfun_glm_prob(final_mods[[3]], newdata = dt_tv)
fwrite(x = dt_p, file = paste0(out_dir, 'pred/', out_base_name, '.csv'))

# compute correlation between model validation predictions as a measure of consistency
cor_p <- cor(x = as.vector(dt_p[tv=="v", rf_prob1]), y=as.vector(dt_p[tv=="v", glm_prob1]), method = 'pearson')
cor_s <- cor(x = as.vector(dt_p[tv=="v", rf_prob1]), y=as.vector(dt_p[tv=="v", glm_prob1]), method = 'spearman')

# save evaluation metrics, including correlation between model predictions
dt_e <- final_mods[[1]]
dt_e$rfglm_pcor <- cor_p
dt_e$rfglm_scor <- cor_s
fwrite(x = dt_e, file = paste0(out_dir, 'eval/', out_base_name, '.csv'))

# save the GLM model coefficients
dt_gcoefs <- as.data.table(summary(final_mods[[3]])$coefficients, keep.rownames = TRUE)
colnames(dt_gcoefs) <- c("pred_var", colnames(dt_gcoefs)[2:5])
fwrite(x = dt_gcoefs, file = paste0(out_dir, 'glm_coefs/', out_base_name_simp, '.csv'))


# get importance and save
# traditional RF importance (using training data)
dt_rfimp <- as.data.table(randomForest::importance(x = final_mods[[2]])[,3:4], keep.rownames = TRUE)
colnames(dt_rfimp) <- c("pred_var", "rfimp_mndacc", "rfimp_mndgin")

# AUC permutation importance using validation data
# make sure pa is a factor in the validation table
dt_v$pa <- as.factor(dt_v$pa)
# RF AUC perm imp
dt_rfimp_pauc <- as.data.table(vi(final_mods[[2]],
                                  method = "permute",
                                  train = dt_v,
                                  target = "pa",
                                  feature_names = pred_var_names_fs,
                                  metric = "roc_auc",
                                  event_level = "second",
                                  pred_wrapper = pfun_rf_prob,
                                  nsim = 100,
                                  type = 'difference'))
colnames(dt_rfimp_pauc) <- c("pred_var", "rfimp_mndauc", "rfimp_sddauc")


# GLM AUC perm imp
dt_glmimp_pauc <- as.data.table(vi(final_mods[[3]],
                                   method = "permute",
                                   train = dt_v,
                                   target = "pa",
                                   feature_names = pred_var_names_fs,
                                   metric = "roc_auc",
                                   event_level = "second",
                                   pred_wrapper = pfun_glm_prob,
                                   nsim = 100,
                                   type = 'difference'))
colnames(dt_glmimp_pauc) <- c("pred_var", "glmimp_mndauc", "glmimp_sddauc")

# join AUC perm tables
dt_imp_pauc <- dt_glmimp_pauc[dt_rfimp_pauc, on = "pred_var"]

# join AUC perm table to RF imp
dt_imp <- dt_imp_pauc[dt_rfimp, on = "pred_var"]

# run Boruta,a good option if the features aren't correlated
# dont exclude anything
keep_cols_b <- c("pa", pred_var_names_fs)
dt_t_b <- dt_t[ , ..keep_cols_b]
bor <- Boruta::Boruta(pa ~ ., dt_t_b, maxRuns = 100)
dt_b <- as.data.table(Boruta::attStats(bor), keep.rownames = TRUE)
colnames(dt_b) <- c("pred_var", colnames(dt_b)[2:6], "boruta_decision")
dt_b <- dt_b[,c(1,2,7)]
colnames(dt_b) <- c("pred_var", "borimp_mnzdacc", "bor_dec")

# join Boruta to other importance metrics
dt_imp <- dt_b[dt_imp, on = "pred_var"]

# join glm coef estimate
dt_gcoefs_e <- dt_gcoefs[,c(1,2)]
colnames(dt_gcoefs_e) <- c("pred_var", "glmimp_coefest")
dt_imp <- dt_gcoefs_e[dt_imp, on = "pred_var"]

dt_imp <- dt_imp[order(rfimp_mndauc, decreasing = TRUE)]
fwrite(x = dt_imp, file = paste0(out_dir, 'imp/', out_base_name, '.csv'))


# PDPs
cat("Calculating partial dependence...\n")
pdp_list <- list()
for(v in 1:length(pred_var_names_fs)){
  rf_pdp_v <- pdp::partial(object = final_mods[[2]], pred.var = pred_var_names_fs[v], which.class="1", type = "classification",
                           train = final_mods[[4]], prob = TRUE, grid.resolution = 101, plot = FALSE)
  
  glm_pdp_v <- pdp::partial(object = final_mods[[3]], pred.var = pred_var_names_fs[v], type = "classification",
                            train = final_mods[[4]],  grid.resolution = 101, plot = FALSE)
  glm_probs <- glm_pdp_v[,2]
  
  pdp_x_q <- ecdf(rf_pdp_v[,1])(rf_pdp_v[,1])
  dt_pdp_v <- data.table(pred_var = pred_var_names_fs[v], pdp_x = rf_pdp_v[,1], pdp_x_q = pdp_x_q, 
                         rf_pdp_y = rf_pdp_v[,2], glm_pdp_y = glm_probs)
  pdp_list[[v]] <- dt_pdp_v
}
dt_pdp_c <- rbindlist(pdp_list)
fwrite(x = dt_pdp_c, file = paste0(out_dir, 'pdp/', out_base_name, '.csv'))


# Export RF trees
# Prep model
rf_prep = prep.mod(final_mods[[2]], 'randomForest', 'regression', 'base', 1)

# Convert forest
cat("Converting RF trees to text file... \n")
convert.forest(rf_prep, paste0(out_dir, 'rftrees/', out_base_name_simp, '.txt'))
closeAllConnections() 


