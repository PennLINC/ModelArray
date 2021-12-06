"
Author: Bart Larse
Date: 12/6/2021

The most basic way to use this is to do:
visualize_model(my_model, smooth_var = 'Age')

If you want to plot an interaction, you can use int_var = 'sex'' for example

show.data is true or false if you want to plot the data points

There are also options to change the plot labels
plabels = NULL, #plot tile
xlabel=NULL, #xaxis label
ylabel=NULL, # y axis label

You can plot the derivative using derivative_plot=TRUE

UNDER THE HOOD: - Chenying's notes:
- using the fitted gam model to predict the values to plot (middle line, se's boundary), but not to re-fit it (except for derivatives)
"

## Load libraries
library(cowplot)
library(gratia)
library(scales)
library(ggside)
library(ggplot2)   # added by Chenying

font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
theme_replace(axis.text=element_text(colour = "black",size = font_size))
line_size <- 1.5
point_size <- 2
### function to extract derivative, confidence interval, significance, and plot for GAMs ###

get_derivs_and_plot <- function(modobj,smooth_var,low_color=NULL,hi_color=NULL,xlabel=NULL){
  this_font_size = font_size
  if (is.null(low_color)){low_color = "white"}
  if (is.null(hi_color)){hi_color = "grey20"}

  derv<-derivatives(modobj,term=sprintf('s(%s)',smooth_var))
  derv<- derv %>%
    mutate(sig = !(0 >lower & 0 < upper))
  derv$sig_deriv = derv$derivative*derv$sig
  cat(sprintf("\nSig change: %1.2f - %1.2f\n",min(derv$data[derv$sig==T]),max(derv$data[derv$sig==T])))
  d1<- ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = sig_deriv))
  
  # Set the gradient colors
  # if (min(derv$sig_deriv)>=0) {
  #   d1 <- d1 + scale_fill_gradient(low = low_color, high = hi_color,limits = c(0,max(derv$sig_deriv)))
  #   # If you want all plots to have the same scaling, this code can be used instead-- This is desirable if you have a factor-smooth model.
  #   # max_val = 5
  #   # d1 <- d1 +scale_fill_gradient(low = low_color,high = hi_color,limits = c(0,max_val),oob = squish)
  # } else if (min(derv$derivative)<0 & max(derv$derivative)<0) {
  # d1 <- d1 + scale_fill_gradient(low = hi_color, high = low_color,limits = c(min(derv$derivative),0))
  # }else {
  d1 <- d1 +
  scale_fill_gradient2(low = "steelblue", midpoint = 0, mid = "white",
  high = "firebrick",limits = c(min(derv$derivative),max(derv$derivative)))
  # }
  
  if (!is.null(xlabel)) {
    d1 <- d1 + 
      labs(x = xlabel,fill = sprintf("\u0394%s",smooth_var))
  } else{
    d1 <- d1 + 
      labs(x = smooth_var,fill = sprintf("\u0394%s",smooth_var))
  }
  d1 <- d1 + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = this_font_size),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=this_font_size),
          legend.text = element_text(size = this_font_size),
          axis.title = element_text(size = this_font_size),
          legend.key.width = unit(1,"cm"),
          legend.position = "bottom",
          legend.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    guides(fill = guide_colorbar(reverse = F,direction = "horizontal",title.position = "left")) +
    geom_rect(aes(ymin=0,ymax=1,xmin=min(data),xmax=max(data)),color="black",fill="white",alpha = 0)
  return(d1)
}

# Func to visualize GAM model outputs
visualize_model <- function(modobj,
                            smooth_var, 
                            int_var = NULL ,
                            group_var = NULL, 
                            plabels = NULL,
                            xlabel=NULL,
                            ylabel=NULL,
                            ymax = NA,
                            quantile_limits = FALSE,
                            quants = c(.025,.975),
                            check_diagnostics = FALSE,
                            derivative_plot = FALSE,
                            difference_plot = FALSE,
                            show.data=TRUE, 
                            line_color = "black",
                            side_density=FALSE,
                            show.legend=TRUE){
  this_font_size = font_size
  if (any(class(modobj)=="gam")) {
    model <- modobj
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  s<-summary(model)
  df <- model$model
  ## Generate custom line plot
  np <- 10000 #number of predicted values
  
  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  
  if (!is.null(int_var)) {
    # We will produce and interaction plot
    if (!any(grepl(x=as.character(model$formula),pattern = int_var))) {
      warning("int_var not recognized in model formula!")
      return()
    }
    switch (varClasses[int_var],
            "numeric" = {
              q <- quantile(df[,int_var],probs = c(.1,.9)) #pick 10% and 90% to plot
              bigq <- q[[2]]
              smallq <- q[[1]]
              values <- c(bigq,smallq)
              labs <- c(sprintf("high (%1.2f)",bigq),sprintf("low (%1.2f)",smallq))
              
              q <-quantile(rescale(df[,int_var],c(0,1)),probs = c(0,.5,1))
              limit_values <- c(q[[1]],q[[length(q)]])
              midpoint_val <- unname(q[[2]])
              cbar_vals <- unname(q)
              
              theseLabs = rep(values,each = np)
              grad_fill = T
            },
            "factor" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = rep(values,each = np)
              grad_fill = F
            },
            "ordered" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = ordered(rep(values,each = np),levels = values)
              grad_fill = F
            }
    )
    labPred <- data.frame(init = rep(0,np*length(labs)))
    labPred[,int_var] = theseLabs
    labPred$lab = rep(labs,each = np)
    labPred <- labPred[,names(labPred) !="init"]
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else if (v == int_var) {
        next
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    
    thisPred <- thisPred %>% select(-init)
    thisPred <- do.call("rbind", replicate(length(labs), thisPred, simplify = FALSE))
    
    pred <- cbind(labPred,thisPred)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    if (!is.null(group_var)) {
      pred[,group_var] = NA #these columns have to exist in the dataframe for plotting
    }
    pred[,thisResp] = 1 #these columns have to exist in the dataframe for plotting
    
    low_color = "#91bfdb"
    high_color = "#fc8d59"
    high_line = "#f46d43"
    low_line = "#4575b4"

    if (grad_fill == T) {
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var)) 
      if (show.data==TRUE) {
        p1<-p1+ geom_point(alpha = 0.65,stroke = 0, size = point_size) 
      }
       
      if (!is.null(group_var)) {
        if (show.data==TRUE) {
          cat("adding lines")
          p1<- p1 + geom_line(aes_string(group = group_var),alpha = .5)
        }
      }
      p1 <- p1 +
        scale_color_gradientn(colors = c(low_line,"grey90",high_line), values = cbar_vals,name = "") +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = "lab"),alpha = .3, linetype = 0) +
        scale_fill_manual(values = c(high_color,low_color)) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",group = "lab"),size = line_size) +
        labs(title = plabels)
      if (!is.null(xlabel)) {
        p1<-p1+xlab(xlabel)+ylab(ylabel)
      }
      if (side_density==TRUE) {
        p1 <- p1 + geom_xsidedensity(aes(y=stat(density)),show.legend = FALSE) +
          geom_ysidedensity(aes(x=stat(density)),show.legend = FALSE)
      }
    } else {
      black_color = scales::muted('blue',l=60,c=80)#'blue' #"#1c405eff"
      green_color = scales::muted('red',l=60,c=80)# "#285c27ff"
      
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var,fill=int_var))
      if (show.data==TRUE) {
        print('show data is on')
        p1<- p1 +  
          geom_point(alpha = .35,stroke = 0, size = point_size,show.legend = FALSE)
          # geom_hex(color=NULL)
      } 
      
      if (!is.null(group_var)) {
        if (show.data==TRUE) {
          p1<- p1 + geom_line(aes_string(group = group_var),alpha = .3)
        }
      } 
      col_labs <- c(black_color,green_color)
      names(col_labs)=labs
      if (difference_plot == TRUE) {
        col_labs<-c(col_labs,"sig"="black","non_sig"=NA)
      }
      p1 <- p1 +
        # scale_color_brewer(type = "qual",palette = "Set1",direction = 1) +
        scale_color_manual(values = col_labs) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = int_var),alpha = .3, linetype = 0,show.legend=FALSE) +
        # scale_fill_brewer(type = "qual",palette = "Set1",direction = 1) +
        scale_fill_manual(values = col_labs) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",color = int_var),size = line_size,show.legend = TRUE) 
      p1<-p1+theme(legend.position = c(.2,.75),legend.background = element_blank()) 

      if (!is.na(ymax)){
        p1 <- p1 + ylim(NA,ymax)
      }
      if (quantile_limits==TRUE){
        
        q_vals = quantile(df[,thisResp],probs = quants)
        p1 <- p1 + ylim(q_vals[1],q_vals[2])
      }
      if (side_density==FALSE) {
        # p1 <- p1 + geom_rug(sides = "br",alpha = .1)
      }
      p1 <- p1 +labs(title = plabels)+theme(legend.title = element_blank())
      if (!is.null(xlabel)) {
        p1<-p1+xlab(xlabel)+ylab(ylabel)
      }
      if (side_density==TRUE) {
        p1 <- p1 + geom_xsidedensity(aes(y=stat(density)),alpha=.25,show.legend = FALSE) +
          geom_ysidedensity(aes(x=stat(density)),alpha=.25,show.legend = FALSE)
      }
      if (show.legend==FALSE) {
        p1<-p1+theme(legend.position = "none")
      } else {
        #p1<-p1+theme(legend.position = c(1,1))
      }
      #ggExtra::ggMarginal(x = smooth_var,y = thisResp,type = "density")
      
    }
  } else {
    
    # No interaction variable, just produce a single line plot
    int_var = "" # This is a lazy solution to making the existing code workable with no int_var.
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    pred <- thisPred %>% select(-init)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    pred[,group_var] = NA
    pred[,thisResp] = 1
    
    # df <- df %>%
    #   gratia::add_partial_residuals(model)
    # df$partial_resids <- unlist(df[,grep(x=names(df),pattern = "s(",fixed = T)])
    
    p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp))
    if (show.data==TRUE) {
      p1<- p1 +  
        geom_point(alpha = .3,stroke = 0, size = point_size,color = line_color)
        # geom_hex(show.legend = TRUE) + scale_fill_gradient(low="white",high=line_color,limits = c(1, 9), oob = scales::squish)
    } 
    if (!is.null(group_var)) {
      if (show.data==TRUE) {
        cat("adding lines")
        p1<- p1 + geom_line(aes_string(group = group_var),alpha = .5)
      }
    }
    p1 <- p1 + geom_ribbon(data = pred,aes_string(x = smooth_var ,y=thisResp, ymin = "selo",ymax = "sehi"),fill = line_color, alpha = .3, linetype = 0) +
      geom_line(data = pred,aes_string(x = smooth_var, y = "fit"),size = line_size,color=line_color)
    if (side_density==FALSE) {
      # p1 <- p1+geom_rug()
    }
    if (side_density==TRUE) {
      p1 <- p1 + geom_xsidedensity(aes(y=stat(density)),alpha=.25,fill=line_color,show.legend = FALSE) +
        geom_ysidedensity(aes(x=stat(density)),alpha=.25,fill=line_color,show.legend = FALSE)
    }
    p1 <- p1 + labs(title = plabels)
    if (!is.na(ymax)){
      p1 <- p1 + ylim(NA,ymax)
    }
    if (quantile_limits==TRUE){
      
      q_vals = quantile(df[,thisResp],probs = quants)
      p1 <- p1 + ylim(q_vals[1],q_vals[2])
    }
    if (!is.null(xlabel)) {
      p1<-p1+xlab(xlabel)+ylab(ylabel)
    }
  }
  
    if (difference_plot == TRUE) {
    ## First attempt with gratia is not quite right since it does not incorporate the mean differences between groups, only the centered smooths.

    # diffs <-gratia::difference_smooths(model = model, smooth = sprintf("s(%s)",smooth_var),n=1000)
    # diffs<-diffs %>% mutate(significant = factor(!between(x = 0,lower = lower,upper = upper),levels = c("TRUE","FALSE"),labels = c("sig","non_sig")))
    # diffs$yvalue <- min(pred$selo)*.5
    # p1 <- p1 + geom_tile(data=diffs,aes_string(x=smooth_var,fill="significant",y="yvalue"),
    #                      show.legend = FALSE,
    #                      color=NA,
    #                      height=(pred$sehi[[1]]-pred$selo[[1]])/5)+
    #   guides(fill="none")

    p_data <- mgcv::plot.gam(model)
    plot_df <- data.frame(x=p_data[[2]]$x,se=p_data[[2]]$se,fit=p_data[[2]]$fit)
    plot_df <- plot_df %>% mutate(se_upper=fit+se,se_lower=fit-se)%>%
      mutate(significant = factor(!between(x=0,lower=se_lower,upper = se_upper),levels = c("TRUE","FALSE"),labels = c("sig","non_sig")))
    names(plot_df)[names(plot_df)=="x"] = smooth_var
    plot_df$yvalue <- min(pred$selo)*.5

    names(plot_df)[grep(names(plot_df),pattern = 'x')]=smooth_var

    p1 <- p1 + geom_tile(data=plot_df,aes_string(x=smooth_var,fill="significant",y="yvalue"),
                         show.legend = FALSE,
                         color=NA,
                         height=(pred$sehi[[1]]-pred$selo[[1]])/5)+
      guides(fill="none")

    }
  if (derivative_plot == T) {
    # We will add a bar that shows where the derivative is significant.
    # First make some adjustments to the line plot.
    p1<- p1+theme(text = element_text(size=this_font_size),
                  axis.text = element_text(size = this_font_size),
                  axis.title.y = element_text(size = this_font_size),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.text = element_text(size = this_font_size),
                  legend.title = element_text(size = this_font_size),
                  axis.title = element_text(size = this_font_size),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  plot.margin = unit(c(.2, .2, 0, .2), "cm")) #Top, left,Bottom, right
    scatter <- list(p1)

    # Now add the plots using the derivative plotting function
    if (any(grepl(x = row.names(s$s.table),pattern =  ":") & grepl(x=row.names(s$s.table),pattern = int_var))) {
      # Factor levels separately if there is an interaction in the model.
      f<-formula(model) # current formula
      fterms <- terms(f)
      fac <- attr(fterms, "factors")
      idx <- which(as.logical(colSums(fac[grep(x=row.names(fac),pattern = int_var),])))
      new_terms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
      new_formula <- formula(new_terms) # Formula without any interaction terms in the model.
      
      #add derivative gradients for each level of the factor
      num_levels <- length(levels(df[,int_var]))
      level_colors <- suppressWarnings(RColorBrewer::brewer.pal(num_levels,"Set1")) #use the same palette as the line plot
      plotlist = vector(mode = "list",length = num_levels+1) # we will be building a list of plots
      plotlist[1] = scatter # first the scatter plot
      level_colors=c(black_color,green_color)
      for (fcount in 1:num_levels) {
        this_level <- levels(df[,int_var])[fcount]
        df$subset <- df[,int_var] == this_level
        this_mod <- gam(formula = new_formula,data = df,subset = subset)
        # this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        if(!is.null(xlabel)){
          this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount],xlabel = xlabel)
        } else{
          this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        }
        
        if (fcount != num_levels & fcount != 1){
          # get rid of redundant junk
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          this_d$theme$legend.text = element_blank()
          this_d$theme$legend.position="none"
        }
        if (fcount == 1) {
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.text = element_blank()
          this_d$theme$legend.position="none"
        }
        if (fcount == num_levels) {
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          legend <- get_legend(this_d)
        }
        this_d$labels$fill=NULL
        plotlist[fcount+1] = list(this_d+theme(legend.position = "none"))
      }
      
      pg<-cowplot::plot_grid(rel_heights = c(16*num_levels,rep(num_levels,num_levels-1),2*num_levels),plotlist = plotlist,align = "v",axis = "lr",ncol = 1)
      final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.15),ncol = 1)
      # final_plot <- pg
      # print(final_plot)
    } else {
      # No need to split
      if (!is.null(xlabel)) {
        d1 <- get_derivs_and_plot(modobj = modobj,smooth_var = smooth_var,xlabel = xlabel)
      } else{
        d1 <- get_derivs_and_plot(modobj = modobj,smooth_var = smooth_var)
      }
      scatter <- list(p1)
      bar <- list(d1+theme(legend.position = "none"))
      legend <- get_legend(d1+theme(legend.position = "bottom",legend.direction = "horizontal"))
      allplots <- c(scatter,bar)
      pg<-cowplot::plot_grid(rel_heights = c(1,.35),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
      final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.15),ncol = 1)
      # print(final_plot)
    }
    
  }    else {
    # No derivative plot
    p1<- p1+theme(text = element_text(size=font_size),
                  axis.text = element_text(size = font_size),
                  legend.text = element_text(size = font_size),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank())
    final_plot<-p1
  }
  
  # print(final_plot)
  if (check_diagnostics == T) {
    cp <- check(b,
                a.qq = list(method = "tnorm",
                            a.cipoly = list(fill = "light blue")),
                a.respoi = list(size = 0.5),
                a.hist = list(bins = 10))
    print(cp)
  }
  return(final_plot)
}
