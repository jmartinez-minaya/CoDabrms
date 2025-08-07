posterior_smooths_plot <- function(mod, smooth, resp, x, lab.x = "", lab.y = "", lab.z = "")
{
  marg_effe <- brms::posterior_smooths(mod, smooth = smooth, resp = resp)
  x <- mod$data[,x]
  med <- apply(marg_effe, 2, mean)
  #plot(med, sim2$sim2)
  if(is.null(dim(x)))
  {
    quants <- apply(marg_effe, 2, quantile, probs = c(0.025, 0.975))
    data_plot <- data.frame(x = x, med = med, t(quants))
    res <- ggplot(data = data_plot) +
      geom_ribbon(aes(x = x, ymin = X2.5., ymax = X97.5.), fill = "gray") +
      geom_line(aes(x = x, y = med), col = "blue") +
      theme_bw() +
      xlab(lab.x) +
      ylab(lab.y) +
      theme(legend.position="bottom")
  }else{
    #Creating the grid
    x1s <- seq(0, 1, length = 50)
    x2s <- seq(0, 1, length = 50)
    pr <- data.frame(x_s2 = rep(x1s, 50), x_s3 = rep(x2s, rep(50, 50)))

    pred <- brms::posterior_smooths(mod, smooth = smooth, resp = resp, newdata = pr)
    med <- apply(pred, 2, mean)
    
    data_pred <- data.frame(pr, med = med)
    res <- ggplot(data_pred, aes(x_s2, x_s3, z = med)) +
      geom_raster(aes(fill = med)) +
      geom_contour(colour = "white") +
      scale_fill_viridis(name = "") +
      #scale_fill_viridis(name = lab.z) +
      theme_bw() +
      xlab(lab.x) +
      ylab(lab.y) +
      theme(legend.position="bottom",
            panel.background = element_blank(),
            legend.key.height = unit(0.4, "cm"),
            legend.key.width = unit(2, "cm")
            #legend.title = element_text(size = 15),
            #legend.text = element_text(size = 13)
            )
      #zlab(lab.z)
    res <- list(plot_res = res, range_pred = range(med))
  }
  res
}




