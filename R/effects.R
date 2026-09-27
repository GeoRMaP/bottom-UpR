# Hierarchical/random effect helpers ------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.normalize_hierarchical_effects <- function(hierarchical_effects = NULL,
                                            random_effects = NULL) {
  if (is.null(hierarchical_effects) && !is.null(random_effects)) {
    hierarchical_effects <- as.list(random_effects)
  }
  if (is.null(hierarchical_effects) || !length(hierarchical_effects)) return(list())

  lapply(seq_along(hierarchical_effects), function(i) {
    z <- hierarchical_effects[[i]]
    if (is.character(z) && length(z) == 1L) {
      return(list(
        column = z, model = "iid", prediction = "known",
        label = z, internal = paste0(".bottom_re", i)
      ))
    }
    if (!is.list(z)) stop("Each hierarchical effect must be a column name or a list.", call.=FALSE)
    column <- z$column %||% z$name
    if (is.null(column) || length(column) != 1L)
      stop("Hierarchical-effect list entries require 'column' (or 'name').", call.=FALSE)
    list(
      column = as.character(column),
      model = as.character(z$model %||% "iid"),
      prediction = match.arg(z$prediction %||% "known", c("known","zero")),
      label = as.character(z$label %||% column),
      internal = paste0(".bottom_re", i),
      constr = z$constr,
      scale.model = z$scale.model
    )
  })
}

.prepare_hierarchical_effects <- function(data, effects, training_maps = NULL) {
  maps <- vector("list", length(effects))
  if (!length(effects)) return(list(data=data, maps=maps))

  for (i in seq_along(effects)) {
    e <- effects[[i]]
    x <- data[[e$column]]
    if (is.null(training_maps)) {
      lev <- unique(as.character(x[!is.na(x)]))
      map <- stats::setNames(seq_along(lev), lev)
    } else {
      map <- training_maps[[i]]$map
    }
    id <- unname(map[as.character(x)])
    data[[e$internal]] <- as.integer(id)
    maps[[i]] <- list(
      column=e$column, internal=e$internal, label=e$label,
      model=e$model, prediction=e$prediction, map=map
    )
  }
  list(data=data,maps=maps)
}

.effect_formula_terms <- function(effects) {
  if (!length(effects)) return(character())
  vapply(effects, function(e) {
    args <- c(sprintf("f(%s", e$internal), sprintf("model='%s'", e$model))
    if (!is.null(e$constr)) args <- c(args, sprintf("constr=%s", if (isTRUE(e$constr)) "TRUE" else "FALSE"))
    if (!is.null(e$scale.model)) args <- c(args, sprintf("scale.model=%s", if (isTRUE(e$scale.model)) "TRUE" else "FALSE"))
    paste0(paste(args, collapse=", "), ")")
  }, character(1))
}
