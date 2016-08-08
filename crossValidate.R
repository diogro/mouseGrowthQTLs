setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

growth_data = inner_join(growth_phen_std, markers, by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])

#growth_data = inner_join(growth_phen_std, markers, by = "ID")
current_data = growth_data
separateData = function(current_data, per_family = 1){
    current_data = ddply(current_data, .(FAMILY), function(x){
              n_ind = nrow(x)
              x$train = 0
              x$train[sample(1:n_ind, per_family)] = 1
              x }) %>%tbl_df
    return(list(test = filter(current_data, train == 1),
                train = filter(current_data, train == 0)))
}
x = separateData(growth_data)
train_growth_data = x$train
test_growth_data = x$test

current_chrom = 7
fit_SUR_HC = stan(file = './SUR_horseShoe.stan',
                         data = getStanInput(current_chrom, train_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)
fit_SUR_HCp = stan(file = './SUR_horseShoePlus.stan',
                         data = getStanInput(current_chrom, train_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)
fit_SUR_HAL = stan(file = './SUR_HAL.stan',
                         data = getStanInput(current_chrom, train_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)
fit_SUR = stan(file = './SUR.stan',
                         data = getStanInput(current_chrom, train_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)

test_data = getStanInput(current_chrom, test_growth_data, growth_traits)
getMeanSqError = function(model, test_data){
    fit = rstan::extract(model)
    getResidual = function(i){
        residual = test_data$y -
            (fit$w0[i,] + fit$beta_family[i,test_data$family,] + t(fit$w_ad[i,,] %*% t(test_data$ad)) + t(fit$w_dm[i,,] %*% t(test_data$dm)))
        sum(residual^2)
    }
    aaply(1:300, 1, getResidual)
}
getMeanSqError(fit_SUR_HC, test_data)
