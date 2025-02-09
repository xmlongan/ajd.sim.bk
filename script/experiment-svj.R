library(ajd.sim.bk)
library(foreach)
library(doParallel)
source("./script/write_to_csv.R")
#---------------------SVJ-------------------------------------------------------
S = 100; K = 100; v0 = 0.008836; k = 3.99; theta = 0.014; sigma = 0.27
rho = -0.79; r = 0.0319; tau = 5; lmbd = 0.11; mu_b = -0.12; sigma_s = 0.15
true_price = 20.1642 # true option price
#
cl = makeCluster(10)
registerDoParallel(cl)
#
N = 10000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_b, sigma_s, true_price)
write_to_csv(err_dur, "./script/svj-moms8-0010K.csv")
#
N = 40000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_b, sigma_s, true_price)
write_to_csv(err_dur, "./script/svj-moms8-0040K.csv")
#
N = 160000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_b, sigma_s, true_price)
write_to_csv(err_dur, "./script/svj-moms8-0160K.csv")
#
N = 640000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_b, sigma_s, true_price)
write_to_csv(err_dur, "./script/svj-moms8-0640K.csv")
#
N = 2560000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_b, sigma_s, true_price)
write_to_csv(err_dur, "./script/svj-moms8-2560K.csv")
#
stopCluster(cl)
