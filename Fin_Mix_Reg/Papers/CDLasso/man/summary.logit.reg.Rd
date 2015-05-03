\name{summary.logit.reg}
\alias{summary.logit.reg}

\title{
Print all results of Cyclic Coordinate Descent for Logistic Regression
}
\description{
Print full summary of results of Cyclic Coordinate Descent for Logistic Regression.  
}
\usage{
\method{summary}{logit.reg}(object, ...)
}

\arguments{
  \item{object}{
Output of logit.reg. Must be of class "logit.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{summary.logit.reg} produces full output from \link{logit.reg}. For selected output, see \link{print.logit.reg}.
}

\value{
  \item{X}{
  The design matrix.
  }
  \item{Y}{
  The outcome variable for cases.
  }
  \item{cases}{
  The number of cases
  }
  \item{predictors}{
  The number of predictors
  }  
  \item{lambda}{
  The value of penalization parameter \code{lambda} used.
  }
  \item{objective}{
  The value of the objective function
  }
  \item{residual}{
  A vector of length \code{p} listing the residuals
  }
  \item{estimate}{
  The estimate of the coefficients
  }
  \item{nonzeros}{
  The number of "selected" variables included in the model.
  }
  \item{selected}{
  The name of the "selected" variables included in the model.
  }
}


\references{
	Wu, T.T., Chen, Y.F., Hastie, T., Sobel E. and Lange, K. (2009). Genome-wide association analysis by lasso penalized logistic regression. Bioinformatics, Volume 25, No 6, 714-721.

}

\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{summary.logit.reg}}

\code{\link{logit.reg}}
}
\examples{
set.seed(1001)
n=500;p=5000
beta=c(1,1,1,1,1,rep(0,p-5))
x=matrix(rnorm(n*p),p,n)
xb = t(x) \%*\% beta
logity=exp(xb)/(1+exp(xb))
y=rbinom(n=length(logity),prob=logity,size=1)

rownames(x)<-1:nrow(x)
colnames(x)<-1:ncol(x)

#Lasso penalized logistic regression
out<-logit.reg(x,y,lambda=50)

#Re-estimate parameters without penalization
out2<-logit.reg(x[out$selected,],y,lambda=0)
summary(out2)
}