\name{Sensory3Phase.dat}
\docType{data}
\alias{Sensory3Phase.dat}
\alias{Sensory3PhaseShort.dat}
\title{Data for the three-phase sensory evaluation experiment in Brien, C.J. and Payne,
R.W. (1999)}
\description{
The data is from an experiment involved two phases.  In the field phase 
a viticultural
experiment was conducted to investigate the differences between 4
types of trellising and 2 methods of pruning.  The design was a
split-plot design in which the trellis types were assigned to the main
plots using two adjacent Youden squares of 3 rows and 4 columns.  Each
main plot was split into two subplots (or halfplots) and the methods
of pruning assigned at random independently to the two halfplots in
each main plot.  The produce of each halfplot was made into a wine so
that there were 24 wines altogether.

The second phase was an evaluation phase in which the produce from the
halplots was evaluated by 6 judges all of whom took part in 24
sittings.  In the first 12 sittings the judges evaluated the wines
made from the halfplots of one square; the final 12 sittings were to
evaluate the wines from the other square.  At each sitting, each judge
assessed two glasses of wine from each of the halplots of one of the
main plots.  The main plots allocated to the judges at each sitting
were determined as follows.  For the allocation of rows, each occasion
was subdivided into 3 intervals of 4 consecutive sittings.  During
each interval, each judge examined plots from one particular row,
these being determined using two 3x3 Latin squares for each occasion,
one for judges 1-3 and the other for judges 4-6.  At each sitting
judges 1-3 examined wines from one particular column and judges 4-6
examined wines from another column.  The columns were randomized to
the 2 sets of judges x 3 intervals x 4 sittings using duplicates of a
balanced incomplete block design for v=4 and k=2 that were latinized.
This balanced incomplete block design consists of three sets of 2
blocks, each set containing the 4 "treatments".  For each interval, a
different set of 2 blocks was taken and each block assigned to two
sittings, but with the columns within the block placed in reverse
order in one sitting compared to the other sitting.  Thus, in each
interval, a judge would evaluate a wine from each of the 4 columns.

The \code{data.frame} contains the following factors, in  the order give: 
Occasion, Judges, Interval, Sittings, Position, Squares, Rows, Columns, 
Halfplot, Trellis, Method. They are followed by the simulated response 
variable Score.

The scores are ordered so that the factors Occasion, Judges, Interval,
Sittings and Position are in standard order; the remaining factors are
in randomized order.

See also the vignette accessed via \code{vignette("DesignNotes", package="dae")}.
}
\usage{data(Sensory3Phase.dat)
data(Sensory3PhaseShort.dat)}
\format{A data.frame containing 576 observations of 12 variables. 
        There are two versions, one with shorter factor names than the other.
}
\references{Brien, C.J. and Payne, R.W. (1999) Tiers, structure formulae 
and the analysis of complicated experiments. \emph{The Statistician}, \bold{48}, 41-52.
}
\keyword{datasets}