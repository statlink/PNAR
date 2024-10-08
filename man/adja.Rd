\name{adja}
\alias{adja}
\title{
Generation of a network from the Stochastic Block Model
}

\description{
This function generates a network from the Stochastic Block Model with \eqn{K} blocks.
}

\usage{
adja(N, K, alpha, directed = FALSE)
}

\arguments{
\item{N}{
The number of nodes on the network.
}
\item{K}{
The number of blocks. Each block has dimension \eqn{N/K}. \eqn{K}
should be chosen such that \eqn{N} is divisible by \eqn{K}.
}
\item{alpha}{
The network density. A value in \eqn{[0,1]} defining the frequency of connections
in the network.
}
\item{directed}{
Logical scalar, whether to generate a directed network or not. If TRUE a directed
network is generated.
}
}

\details{
For each pair of nodes it performs a Bernoulli trial with values 1 "draw an edge", 0 "otherwise".
The probabilities of these trials are bigger if the two nodes are in the same block, lower otherwise, and they are specified based on the number of nodes on the network \eqn{N} and network density \eqn{alpha}:
Probability to draw an edge for a pair of nodes in the same block: \eqn{\alpha*N^{-0.3}}.
Probability to draw an edge for a pair of nodes in different blocks: \eqn{\alpha*N^{-1}}.
}


\value{
A row-normalized non-negative matrix describing the network. The main diagonal
entries of the matrix are zeros, all the other entries are non-negative and the
sum of elements over the rows equals one.
}

\references{
Faust, K. and S. Wasserman (1992). Blockmodels: Interpretation and evaluation.
\emph{Social Networks}, 14, 5--61.
}

\seealso{
\code{\link{adja_gnp}}
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\examples{
W <- adja(N = 20, K = 5, alpha = 0.1)
}


