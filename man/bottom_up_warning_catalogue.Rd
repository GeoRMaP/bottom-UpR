\name{bottom_up_warning_catalogue}
\alias{bottom_up_warning_catalogue}
\title{List bottom-UpR warning messages and recommended actions}
\usage{bottom_up_warning_catalogue()}
\value{A data frame with warning identifier, category, message pattern, trigger, and suggested action.}
\description{Returns the complete catalogue of package-level warning messages or message patterns emitted by the bottom-UpR workflow. Dynamic values are shown in braces. Warnings originating internally from INLA are not reproduced, although bottom-UpR warnings about INLA diagnostic states are included.}
