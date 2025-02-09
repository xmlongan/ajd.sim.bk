library(ajd.sim.bk)
library(foreach)
library(doParallel)
source("./script/write_to_csv.R")
#---------------------SVCJ------------------------------------------------------
S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1; sigma_s = 0.0001
mu_v = 0.05; rhoJ = -0.38
true_price = 6.8619 # true option price
#
# mu = r - lmbd * mu_b
# mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
h = 1
#
cl = makeCluster(12)
registerDoParallel(cl)
#
N = 10000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho,
  lmbd, mu_b, sigma_s, mu_v, rho_J, true_price)
write_to_csv(err_dur, "./script/svcj-BK-0010K.csv")
#
N = 40000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho,
  lmbd, mu_b, sigma_s, mu_v, rho_J, true_price)
write_to_csv(err_dur, "./script/svcj-BK-0040K.csv")
#
N = 160000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho,
  lmbd, mu_b, sigma_s, mu_v, rho_J, true_price)
write_to_csv(err_dur, "./script/svcj-BK-0160K.csv")
#
N = 640000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho,
  lmbd, mu_b, sigma_s, mu_v, rho_J, true_price)
write_to_csv(err_dur, "./script/svcj-BK-0640K.csv")
#
N = 2560000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.bk::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho,
  lmbd, mu_b, sigma_s, mu_v, rho_J, true_price)
write_to_csv(err_dur, "./script/svcj-BK-2560K.csv")
#
stopCluster(cl)
