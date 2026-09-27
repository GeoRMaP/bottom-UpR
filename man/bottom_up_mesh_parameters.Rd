\name{bottom_up_mesh_parameters}
\alias{bottom_up_mesh_parameters}
\title{Derive data-adaptive SPDE mesh parameters}
\usage{
bottom_up_mesh_parameters(data, coords = c("x", "y"),
  inner_edge_fraction = 0.08, outer_edge_fraction = 0.16,
  inner_offset_fraction = 0.05, outer_offset_fraction = 0.10,
  cutoff_fraction = 0.005, min_cutoff = 0)
}
\arguments{
\item{data}{Data frame containing projected coordinates.}
\item{coords}{Two projected coordinate columns.}
\item{inner_edge_fraction}{Fraction of the study-area diagonal for the inner maximum edge.}
\item{outer_edge_fraction}{Fraction of the study-area diagonal for the outer maximum edge.}
\item{inner_offset_fraction}{Fraction of the study-area diagonal for the inner extension.}
\item{outer_offset_fraction}{Fraction of the study-area diagonal for the outer extension.}
\item{cutoff_fraction}{Fraction of the study-area diagonal used as cutoff.}
\item{min_cutoff}{Optional lower bound for cutoff in coordinate units.}
}
\value{A list of derived mesh distances and spatial-footprint information.}
\description{Scales mesh distances to the location and geographic extent of the supplied projected coordinates.}
