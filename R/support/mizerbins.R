# Mizer Bin Structure code from sizeSpectra::eightMethods


# Default values
# function (oneYear = 1980, 
#           dataForYear = dplyr::filter(data, Year == oneYear), 
#           figName = "eightMethods") 

x         = rep(dataForYear$bodyMass, dataForYear$Number)
log.x     = log(x)
sum.log.x = sum(log.x)
xmin      = min(x)
xmax      = max(x)
num.bins  = 8




# Method 1 - set table structure here:
hLlin.list    = Llin.method(x, num.bins = num.bins)
hLlin.b       = hLlin.list$slope
hLlin.confMin = hLlin.list$confVals[1]
hLlin.confMax = hLlin.list$confVals[2]
eightMethodsRes = data.frame(Year = oneYear, 
                             Method = "Llin", 
                             b = hLlin.b, 
                             confMin = hLlin.confMin, 
                             confMax = hLlin.confMax, 
                             row.names = NULL)



# Skip to Mizer Code

# Same bin number k = 8
hLBmiz.num.bins = num.bins

# Solve for beta (value for equal bin widths) using nlm
beta = nlm(LBmizbinsFun, 2, xmin = xmin, xmax = xmax, k = hLBmiz.num.bins)$est

# Then get bins
hLBmiz.bins = c( beta^(0:(hLBmiz.num.bins - 1)) * min(x), max(x))

# Plotting the histogram
hLBmiz = hist(x, breaks = hLBmiz.bins, plot = FALSE)

# getting minimum of bins (left hand value of bin)
hLBmiz.log.min.of.bins = log(hLBmiz.bins[-length(hLBmiz.bins)])

# getting log of counts
hLBmiz.log.counts = log(hLBmiz$counts)

# inf to NA
hLBmiz.log.counts[is.infinite(hLBmiz.log.counts)] = NA

# calculate slope
hLBmiz.lm = lm(hLBmiz.log.counts ~ hLBmiz.log.min.of.bins, 
               na.action = na.omit)

# pull out slope coefficient
hLBmiz.b = hLBmiz.lm$coeff[2] - 1

# 95% conf
hLBmiz.conf = confint(hLBmiz.lm, "hLBmiz.log.min.of.bins", 0.95) - 1

# Add to table for 8 methods
eightMethodsRes = rbind(eightMethodsRes, 
                        data.frame(Year = oneYear, 
                                   Method = as.factor("LBmiz"), 
                                   b = hLBmiz.b, 
                                   confMin = hLBmiz.conf[1], 
                                   confMax = hLBmiz.conf[2], 
                                   row.names = NULL))