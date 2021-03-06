\docType{data}
\name{tbc}
\alias{tbc}
\alias{tbc.target}
\alias{terneuzen}
\title{Terneuzen birth cohort}
\format{\code{tbs} is a data frame with 3951 rows and 11 columns:
\describe{
\item{id}{Person number}
\item{occ}{Occasion number}
\item{nocc}{Number of occasions}
\item{first}{Is this the first record for this person? (TRUE/FALSE)}
\item{typ}{Type of data (all observed)}
\item{age}{Age (years)}
\item{sex}{Sex 1=M, 2=F}
\item{hgt.z}{Height Z-score}
\item{wgt.z}{Weight Z-score}
\item{bmi.z}{BMI Z-score}
\item{ao}{Adult overweight (0=no, 1=yes)}
}

\code{tbc.target} is a data frame with 2612 rows and 3 columns:
\describe{
\item{id}{Person number}
\item{ao}{Adult overweight (0=no, 1=yes)}
\item{bmi.z.jv}{BMI Z-score as young adult (18-29 years)}
}}
\source{
  De Kroon, M. L. A., Renders, C. M., Kuipers, E. C., van
  Wouwe, J. P., van Buuren, S., de Jonge, G. A., Hirasing,
  R. A. (2008). Identifying metabolic syndrome without
  blood tests in young adults - The Terneuzen birth cohort.
  \emph{European Journal of Public Health}, \emph{18}(6),
  656-660.

  De Kroon, M. L. A., Renders, C. M., Van Wouwe, J. P., Van
  Buuren, S., Hirasing, R. A. (2010).  The Terneuzen birth
  cohort: BMI changes between 2 and 6 years correlate
  strongest with adult overweight.  \emph{PLoS ONE},
  \emph{5}(2), e9155.

  De Kroon, M. L. A. (2011).  \emph{The Terneuzen Birth
  Cohort. Detection and Prevention of Overweight and
  Cardiometabolic Risk from Infancy Onward.} Disseration,
  Vrije Universiteit, Amsterdam.
  \url{http://dare.ubvu.vu.nl/handle/1871/23806}

  Van Buuren, S. (2012). \emph{Flexible Imputation of
  Missing Data.} Boca Raton, FL: Chapman \& Hall/CRC Press.
}
\description{
  Data of subset of the Terneuzen Birth Cohort data on
  child growth.
}
\details{
  This \code{tbc} data set is a random subset of persons
  from a much larger collection of data from the Terneuzen
  Birth Cohort. The total cohort comprises of 2604 unique
  persons, whereas the subset in \code{tbc} covers 306
  persons. The \code{tbc.target} is an auxiliary data set
  containing two outcomes at adult age. For more details,
  see De Kroon et al (2008, 2010, 2011).  The imputation
  methodology is explained in Chapter 9 of Van Buuren
  (2012).
}
\examples{
data <- tbc
md.pattern(data)
}
\keyword{datasets}

