setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = chains)
getSqError = function(fit, test_data)
{
    getResidual = function(i){
        residual = test_data$y -
            (fit$w0[i,] + fit$beta_family[i,test_data$family,] +
             t(fit$w_ad[i,,] %*% t(test_data$ad)) +
             t(fit$w_dm[i,,] %*% t(test_data$dm)))
        sum(residual^2)
    }
    n = nrow(fit[[1]])
    aaply(1:n, 1, getResidual)
}
getMeanSqError = function(fit, test_data)
{
        residual = test_data$y -
            (colMeans(fit$w0) + colMeans(fit$beta_family)[test_data$family,] +
             t((colMeans(fit$w_ad)) %*% t(test_data$ad)) +
             t((colMeans(fit$w_dm)) %*% t(test_data$dm)))
        sum(residual^2)
}
getMeanSqErrorShrink = function(fit, test_data, cut_off = 0.5)
{
        ad_mask = colMeans(fit$shrink_ad)
        ad_mask[] = as.numeric(colMeans(fit$shrink_ad)[] > cut_off)
        dm_mask = colMeans(fit$shrink_dm)
        dm_mask[] = as.numeric(colMeans(fit$shrink_dm)[] > cut_off)
        residual = test_data$y -
            (colMeans(fit$w0) + colMeans(fit$beta_family)[test_data$family,] +
             t((ad_mask * colMeans(fit$w_ad)) %*% t(test_data$ad)) +
             t((dm_mask * colMeans(fit$w_dm)) %*% t(test_data$dm)))
        sum(residual^2)
}
separateData = function(current_data, per_family = 1)
{
    current_data = ddply(current_data, .(FAMILY), function(x){
              n_ind = nrow(x)
              x$train = 0
              x$train[sample(1:n_ind, per_family)] = 1
              x }) %>%tbl_df
    return(list(test = filter(current_data, train == 1),
                train = filter(current_data, train == 0)))
}

necropsy_data = inner_join(necropsy_phen_std, markers, by = "ID")
necropsy_data[necropsy_traits] = scale(necropsy_data[necropsy_traits])

crossValidate = function(current_data, trait_vector, iter, warmup){
    x = separateData(current_data)
    train_data = x$train
    test_data = getStanInputFullGenome(x$test, trait_vector)

    model_list = list(SUR = "./SUR.stan",
                      HC = "./SUR_horseShoe.stan",
                      HCp = "./SUR_horseShoePlus.stan",
                      HAL = "./SUR_HAL.stan")
    fits = llply(model_list,
            function(model)
            rstan::extract(stan(file = model,
                data = getStanInputFullGenome(train_data, trait_vector),
                chain=chains,
                iter = iter,
                warmup = warmup)), .parallel = TRUE)
    list(mse = ldply(fits, getMeanSqError, test_data),
         fits = fits,
         separated_data = x)
}
library(doMC)
registerDoMC(10)
chains = 1
crossValidate(necropsy_data, necropsy_traits, iter = 200, warmup = 100)
CV_msqe = llply(1:10, function(x) crossValidate(necropsy_data, necropsy_traits,
                                                iter = 200, warmup = 100),
                .parallel = TRUE)
write_rds(CV_msqe, paste0(Rdatas_folder, "cross_validation_necropsy.rds"))
llply(CV_msqe, '[[', 1) %>% laply('[[',"V1") %>% aaply(1, function(x) x - min(x))


fit = CV_msqe[[1]]$fits$HCp

CV_shrink = laply(CV_msqe, function(part) aaply(seq(0.1, 0.4, 0.01), 1, function(x) getMeanSqErrorShrink(part$fits$HCp,
               getStanInputFullGenome(part$separated_data$test,
               necropsy_traits), x)))
CV_rank = aaply(CV_shrink, 1, rank)
colnames(CV_shrink) = colnames(CV_rank) = seq(0.1, 0.4, 0.01)
