# Utility functions

# Split dataset into training and validation set, with indices returned for 
# cross-validation
split.data <- function(data.size, cv.folds, types) {
  
	type.names = sort(unique(types))
	type.size = floor(table(types)/cv.folds)[type.names]
  
  if(any(type.size == 0)) 
  	cat('Not enough observations in one group for cross-validation!')
  
  index.set = 1:data.size
  
  fold.size = floor(data.size/cv.folds)
  
  train = list()
  valid = list()
  
  for(fold.i in 1:cv.folds) {
  	valid[[fold.i]] = vector()
  	  
  	for(type.i in 1:length(type.names)) {
  		valid[[fold.i]] = c(valid[[fold.i]], 
  												sample(index.set[types[index.set]==type.names[type.i]],
  															 type.size[type.i]))
  	}
  	
    train[[fold.i]] = which(!((1:data.size) %in% valid[[fold.i]]))
    
    index.set = index.set[!(index.set %in% valid[[fold.i]])]
  }
  
  return(list(valid=valid, train=train))
}

# Scale dataset for test or validation
scaleData <- function(X, Y, groups, group.names, x.scale, y.scale, 
                      exclude=rep(FALSE, dim(X)[2]),
                      scale.factors=NULL) {
  
  X.new = X[,!exclude]
  Y.new = Y
  
  p.out = dim(X)[2]
  k = length(group.names)
  
  if(is.null(scale.factors)) {
    x.means = matrix(NA, p.out, k+1)
    x.sds = matrix(NA, p.out, k+1)
    y.mean = rep(NA, k+1)
    y.sd = rep(NA, k+1)
  } else {
    x.means = scale.factors$x.means
    x.sds = scale.factors$x.sds
    y.mean = scale.factors$y.mean
    y.sd = scale.factors$y.sd
  }
  
  # Scaling data group-wise
  for(group.i in 1:k) {
    group.inds = groups == group.names[group.i]
    
    if(sum(group.inds) > 1) {
      # Scale both x and y by y sd for each group
      
      if(is.null(scale.factors)) {
        
        if(sd(Y[group.inds])==0) { # Special case where all cell lines resistant
          Y.temp = scale(Y[group.inds], scale=FALSE) 
          y.sd[group.i] = 1 # to avoid problems when scaling validation data
        } else{
          Y.temp = scale(Y[group.inds]) 
          y.sd[group.i] = attr(Y.temp, 'scaled:scale')
        }
        
        y.mean[group.i] = attr(Y.temp, 'scaled:center')
        
      
        X.temp = scale(X[group.inds,!exclude]) 
        x.means[!exclude,group.i] = attr(X.temp, 'scaled:center')
        x.sds[!exclude,group.i] = attr(X.temp, 'scaled:scale')
        
        if(any(x.sds[!exclude,group.i]==0)) {
          warning(paste('SD 0 in group', group.i))
          X.temp[,x.sds[!exclude,group.i]==0] = 0 # Avoid NAs for identical values
          x.sds[x.sds[!exclude,group.i]==0,group.i] = 1 # Avoid division by zero when using pre-calculated SDs
        }
      } else {
        Y.temp = scale(Y[group.inds], y.mean[group.i], y.sd[group.i])
        X.temp = scale(X[group.inds,!exclude], x.means[!exclude,group.i], 
                       x.sds[!exclude,group.i]) 
      }
      
      if(x.scale == 'local') {
        X.new[group.inds,] = X.temp
      }
      
      if(y.scale == 'local') {
        Y.new[group.inds] = Y.temp
      }
    }
  }
  
  # Scaling data globally
  if(is.null(scale.factors)) {
    
    X.temp = scale(X[,!exclude])
    x.means[!exclude,k+1] = attr(X.temp, 'scaled:center')
    x.sds[!exclude,k+1] = attr(X.temp, 'scaled:scale')
    
    if(any(x.sds[!exclude,k+1]==0)) {
      warning('SD 0 overall.')
      X.temp[,x.sds[!exclude,k+1]==0] = 0 # Avoid NAs for identical values
      x.sds[x.sds[!exclude,k+1]==0,k+1] = 1 # Avoid division by zero when using pre-calculated SDs
    }
    
    Y.temp = scale(Y) 
    y.mean[k+1] = attr(Y.temp, 'scaled:center')
    y.sd[k+1] = attr(Y.temp, 'scaled:scale')
  } else {
    Y.temp = scale(Y, y.mean[k+1], y.sd[k+1])
    X.temp = scale(X[,!exclude], x.means[!exclude,k+1], 
                   x.sds[!exclude,k+1]) 
  }
  
  if(x.scale == 'global') {
    X.new = X.temp
  }
  
  if(y.scale == 'global') {
    Y.new = Y.temp
  }
  
  return(list(X=cbind(X.new, X[,exclude]), Y=Y.new, 
              scale.factors=list(x.means=x.means, x.sds=x.sds,
                                                   y.mean=y.mean, y.sd=y.sd)))
  
}

# Rescaling the response from local to unscaled
rescaleResponseUnscaled <- function(Y, group.i, scale.factors) {
  
  y.mean = scale.factors$y.mean
  y.sd = scale.factors$y.sd
  
  Y.unscaled = Y*y.sd[group.i]+y.mean[group.i]
  
  return(Y.unscaled)
}
