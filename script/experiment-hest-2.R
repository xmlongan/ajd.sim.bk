library(ajd.sim.bk)
library(foreach)
library(doParallel)
source("./script/write_to_csv.R")
#---------------------Heston----------------------------------------------------
# parameter setting 2
S = 100; K = 100; v0 = 0.09; k = 2.00; theta = 0.09; sigma = 1.00
rho = -0.3; r = 0.05; tau = 5
true_price = 34.9998
#
price_hest(1000, S, K, v0, tau, r, k, theta, sigma, rho, true_price, log.inv=T)
#
cl = makeCluster(10)
registerDoParallel(cl)
#
N = 10000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_hest(
  N, S, K, v0, tau, r, k, theta, sigma, rho, true_price)
write_to_csv(err_dur, "./script/heston-s2-BK-0010K.csv")
#
N = 40000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_hest(
  N, S, K, v0, tau, r, k, theta, sigma, rho, true_price)
write_to_csv(err_dur, "./script/heston-s2-BK-0040K.csv")
#
N = 160000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_hest(
  N, S, K, v0, tau, r, k, theta, sigma, rho, true_price)
write_to_csv(err_dur, "./script/heston-s2-BK-0160K.csv")
#
N = 640000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_hest(
  N, S, K, v0, tau, r, k, theta, sigma, rho, true_price)
write_to_csv(err_dur, "./script/heston-s2-BK-0640K.csv")
#
N = 2560000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_hest(
  N, S, K, v0, tau, r, k, theta, sigma, rho, true_price)
write_to_csv(err_dur, "./script/heston-s2-BK-2560K.csv")
#
stopCluster(cl)
