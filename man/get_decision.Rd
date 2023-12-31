% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_decision.R
\name{get_decision}
\alias{get_decision}
\title{Perform A+B+C test for The Rule-Based Design for Agents with Non-Overlapping Dose-Limiting Toxicities}
\usage{
get_decision(d, p.true, p.conditional1, p.conditional2, A, B, C, a.E, a.D, b.E, b.D, c.E)
}
\arguments{
\item{d}{the current dose combination}

\item{p.true}{a list of true DLT rates at all dose combination levels}

\item{p.conditional1}{The conditional probability of observing only DLT1, given that at least one DLT occurs and that DLT3 is absent}

\item{p.conditional2}{The conditional probability of observing only DLT2, given that at least one DLT occurs and that DLT3 is absent}

\item{A}{The sample size of first cohort of A+B+C design}

\item{B}{The sample size of first cohort of A+B+C design}

\item{C}{The sample size of first cohort of A+B+C design}

\item{a.E}{is the maximum number of DLTs needed out of A patients to to conclude the current dose is accepted;}

\item{a.D}{is the minimum number of DLTs needed out of A patients to to conclude the current dose is rejected;}

\item{b.E}{is the maximum number of DLTs needed out of B patients to to conclude the current dose is accepted;}

\item{b.D}{is the minimum number of DLTs needed out of B patients to to conclude the current dose is rejected;}

\item{c.E}{is the minimum number of DLTs needed out of C patients to to conclude the current dose is rejected, otherwise, the current dose will be accepted;}
}
\value{
\code{get_decision()} returns the result as a list,
        including:
        (1) \code{$flag.accepted}: A Boolean value indicating whether the current dose combination should be accepted
        (2) \code{$t.table}: A matrix containing the number of DLT1, DLT2, and DLT3, respectively.
        (3) \code{$n.DLTs}: The number of DLTs happens during the test
        (4) \code{$n.patients}: The number of patients enrolled for the test
}
\description{
Perform A+B+C test for The Rule-Based Design for Agents with Non-Overlapping Dose-Limiting Toxicities
}
\seealso{
Lee BL, Fan SK. A two-dimensional search algorithm for dose-finding trials of two agents. J Biopharm Stat. 2012;22(4):802-18. doi: 10.1080/10543406.2012.676587. PMID: 22651116.
}
