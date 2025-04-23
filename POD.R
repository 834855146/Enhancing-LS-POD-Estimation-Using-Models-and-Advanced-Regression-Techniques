calculate_POD <- function(beta0_mu , beta0_sigma, beta1_mu , beta1_sigma , cov_coef, sigma_hat , var_sigma_hat, th , p=0.9 , q=0.95, two.side=FALSE, t=FALSE){
  if (t) {
    zp = qt(p,9)
  }else{
    zp = qnorm(p)
  }
  
  if (two.side) {
    q = 1 - (1 - q)/2
  }
  u_pod <- (th - beta0_mu)/beta1_mu
  sigma_pod <- sigma_hat/beta1_mu  
  
  v_u <- (beta0_sigma**2 + 2*u_pod*cov_coef + u_pod*beta1_sigma**2)/beta1_mu**2
  cov_u <- (sigma_pod*cov_coef + sigma_pod*u_pod*beta1_sigma**2)/beta1_mu**2
  v_pod <- (var_sigma_hat + sigma_pod**2*beta1_sigma**2)/beta1_mu**2
  
  v_xp <- v_u + 2*zp*cov_u + zp**2*v_pod
  
  ap <- u_pod + zp*sigma_pod
  
  if (two.side) {
    apq_low = ap - qnorm(q)*sqrt(v_xp)
    apq_up = ap + qnorm(q)*sqrt(v_xp)
    return(c(ap,apq_low , apq_up))
  } 
  apq <- ap + qnorm(q)*sqrt(v_xp)
  
  return(c(ap,apq,u_pod))
}

calulate_POD_cruve <- function(beta0_mu , beta0_sigma, beta1_mu , beta1_sigma , cov_coef, sigma_hat , var_sigma_hat, th , p=0.9 , q=0.95, two.side=FALSE , group){
  probs <- seq(0.001, 0.999, by = 0.001)
  POD <- matrix(0, nrow = length(probs)) 
  
  if (two.side) {
    POD_CI_low <- matrix(0, nrow = length(probs)) 
    POD_CI_up <- matrix(0, nrow = length(probs)) 
    
    for (i in 1:length(probs)) {
      res <- calculate_POD(beta0_mu , beta0_sigma, beta1_mu , beta1_sigma , cov_coef, sigma_hat , var_sigma_hat, th,  probs[i], q ,two.side)
      
      POD[i] = res[1]
      POD_CI_low[i] = res[2]
      POD_CI_up[i] = res[3]
    }
    
    result = data.frame(y=c(probs,probs,probs),x=c(POD,POD_CI_low , POD_CI_up),
                        group = group,
                        linetype =  c(factor(rep('ap',length(probs))),
                                      factor(rep('ap/q_low',length(probs))),
                                      factor(rep('ap/q_up',length(probs))))
    )
  }else{
    POD_CI <- matrix(0, nrow = length(probs)) 
    
    for (i in 1:length(probs)) {
      res <- calculate_POD(beta0_mu , beta0_sigma, beta1_mu , beta1_sigma , cov_coef, sigma_hat , var_sigma_hat, th,  probs[i], q ,two.side)
      
      POD[i] = res[1]
      POD_CI[i] = res[2]
    }
    
    result = data.frame(y=c(probs,probs),x=c(POD,POD_CI),
                        group = group,
                        linetype =  c(factor(rep('ap',length(probs))),
                                      factor(rep('ap/q',length(probs))))
    )
  }
  return(result)
}

POD_plt <- function(result , color){
  ggplot(result, aes(x = x, y = y, color = group, group = linetype, linetype = ifelse(grepl("ap/q", linetype), "ap/q", "ap"))) +
    geom_line(linewidth=0.9) +
    xlim(0, max(result$x))+
    labs(x = "a", y = "p") +
    scale_linetype_manual(
      values = c("ap" = "solid", "ap/q" = "dashed")
    ) +
    scale_color_manual(
      values = color
    )+
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 14),   
          legend.key.size = unit(1.5, "cm"))+
    guides(linetype = "none") +
    geom_hline(aes(yintercept=0.9))
}


# function for pod

linear_pod_cal <-function(samples){
  model1 = lm(ahat~a, data = samples)
  
  beta0_mu = model1$coefficients[1]
  beta0_sigma = summary(model1)$coefficients[1,2]
  beta_mu1 = model1$coefficients[2]
  beta_sigma1 = summary(model1)$coefficients[2,2]
  p <- length(coef(model1))
  sigma_hat1 <- summary(model1)$sigma
  var_sigma_hat1 <- sigma_hat1^2 / (2 * (n - p))

  pod = calculate_POD(beta0_mu , beta0_sigma, beta_mu1 , beta_sigma1, 0, sigma_hat1, var_sigma_hat1, th )
  
  return(pod)
}

robust_pod_cal <-function(samples,mod = FALSE){
  if (mod) {
    model2 = rlm(ahat ~ a.mode-1, data = samples, method = "M", maxit = 100,psi=psi.bisquare,c=2)
    
  }else{
    model2 = rlm(ahat ~ a-1, data = samples, method = "M")

  }
  beta_mu2 = model2$coefficients[1]
  beta_sigma2 = summary(model2)$coefficients[1,2] 
  
  p2 <- length(coef(model2))
  sigma_hat2 <- summary(model2)$sigma
  var_sigma_hat2 <- sigma_hat2^2 / (2 * (n - p2))
  
  pod = calculate_POD(0 , 0, beta_mu2 , beta_sigma2 , 0, sigma_hat2 , var_sigma_hat2, th, two.side = FALSE)
  return(pod)
}

physics_pod_cal <-function(samples,b=0.5){
  
  model2 = lm(ahat~a.mode-1, data = samples)
  
  beta_mu2 = model2$coefficients[1]
  beta_sigma2 = summary(model2)$coefficients[1,2]
  
  p2 <- length(coef(model2))
  sigma_hat2 <- summary(model2)$sigma
  var_sigma_hat2 <- sigma_hat2^2 / (2 * (n - p2))
  
  pod = calculate_POD(0 , 0, beta_mu2 , beta_sigma2 , 0, sigma_hat2, var_sigma_hat2, th, two.side = FALSE)^(1/b)
  
  return(pod)
}

weight_pod_cal <-function(samples,mod = FALSE){
  w =  (samples$a)
  if (mod) {
    model_wls <- lm(ahat ~ a.mode-1, data = samples, weights = w)
  }else{
    model_wls <- lm(ahat ~ a-1, data = samples, weights = w)
  }
  
  beta_mu2 = model_wls$coefficients[1]
  beta_sigma2 = summary(model_wls)$coefficients[1,2]
  
  p2 <- length(coef(model_wls))
  sigma_hat2 <- summary(model_wls)$sigma
  var_sigma_hat2 <- sigma_hat2^2 / (2 * (n - p2))
  
  pod = calculate_POD(0 , 0, beta_mu2 , beta_sigma2 , 0, sigma_hat2, var_sigma_hat2, th, two.side = FALSE)
  
  return(pod)
}

bayes_pod_cal <-function(samples,compiled_model,stan_data=stan_data){
  model3 <- sampling(compiled_model, data = stan_data, iter = 2000, chains = 4,refresh = 0)
  sum3 = summary(model3)$summary
  
  beta_mu3 = sum3[1,1]
  beta_sigma3 = sum3[1,3]
  sigma_hat3 <- sum3[2,1]
  var_sigma_hat3 <- sum3[2,3]**2
  
  pod = calculate_POD(0 , 0, beta_mu3 , beta_sigma3 , 0, sigma_hat3 , var_sigma_hat3, th, two.side = FALSE) 
  
  return(pod)
} 

nls_pod_cal <-function(samples , a=0.3 , b = 0.4){
  
  model_nls <- nls(samples$ahat ~ a * samples$a^b,data = samples, start = list(a = a, b = b))
  
  beta_mu_nls = summary(model_nls)$coefficients[1,1]
  beta_sigma_nls = summary(model_nls)$coefficients[1,2]
  sigma_hat_nls <- summary(model_nls)$sigma
  p <- length(coef(model_nls))
  var_sigma_hat_nls <- sigma_hat_nls^2 / (2 * (n - p))
  b = coef(model_nls)[2]
  
  pod = calculate_POD(0 , 0, beta_mu_nls , beta_sigma_nls, 0, sigma_hat_nls, var_sigma_hat_nls, th)^(1/b)
  
  return(list(pod = pod , b = b))
}