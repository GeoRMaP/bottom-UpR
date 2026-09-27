library(shiny)
library(leaflet)

set.seed(2027)

make_comparison_grid <- function() {
  g <- expand.grid(
    lon = seq(8.4, 15.8, length.out = 42),
    lat = seq(2.2, 12.8, length.out = 48)
  )
  n <- nrow(g)

  g$region <- cut(
    g$lat,
    breaks = quantile(g$lat, probs = seq(0, 1, length.out = 6)),
    include.lowest = TRUE,
    labels = paste0("Region ", 1:5)
  )

  g$buildings <- rpois(
    n,
    lambda = pmax(
      0.3,
      exp(
        1.8 +
          0.25 * sin(g$lon * 1.3) +
          0.20 * cos(g$lat * 0.8)
      )
    )
  )

  eta_a <- 0.95 +
    0.18 * sin(g$lon * 1.1) -
    0.14 * cos(g$lat * 0.9)

  eta_b <- 1.00 +
    0.14 * sin(g$lon * 1.0) -
    0.11 * cos(g$lat * 0.8) +
    0.05 * sin((g$lon + g$lat) * 0.7)

  g$model_a <- g$buildings * exp(eta_a)
  g$model_b <- g$buildings * exp(eta_b)
  g$difference <- g$model_b - g$model_a
  g$pct_difference <- 100 * g$difference / pmax(g$model_a, 1e-8)

  g$sd_a <- pmax(0.5, 0.16 * g$model_a + 0.2)
  g$sd_b <- pmax(0.5, 0.18 * g$model_b + 0.2)
  g$cv_a <- g$sd_a / pmax(g$model_a, .Machine$double.eps)
  g$cv_b <- g$sd_b / pmax(g$model_b, .Machine$double.eps)

  zero <- g$buildings == 0
  g[zero, c(
    "model_a", "model_b", "difference", "pct_difference",
    "sd_a", "sd_b", "cv_a", "cv_b"
  )] <- 0

  g
}

demo <- make_comparison_grid()

ui <- fluidPage(
  titlePanel("bottom-UpR: Gridded Model Comparison"),
  sidebarLayout(
    sidebarPanel(
      helpText(
        "Synthetic demo comparing two posterior gridded population surfaces. ",
        "Use this layout to inspect differences between PPB and COUNT models, ",
        "alternative likelihoods, or sensitivity analyses."
      ),
      selectInput(
        "metric",
        "Map variable",
        choices = c(
          "Model A population" = "model_a",
          "Model B population" = "model_b",
          "B - A difference" = "difference",
          "Percent difference" = "pct_difference",
          "Model A CV" = "cv_a",
          "Model B CV" = "cv_b"
        ),
        selected = "difference"
      ),
      selectInput(
        "region",
        "Region",
        choices = c("All", levels(demo$region)),
        selected = "All"
      ),
      sliderInput(
        "radius",
        "Display radius",
        min = 2,
        max = 8,
        value = 4,
        step = 1
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Difference map", leafletOutput("map", height = 650)),
        tabPanel("Model scatter", plotOutput("scatter", height = 550)),
        tabPanel("Regional totals", tableOutput("regional")),
        tabPanel("Diagnostics", plotOutput("diagnostics", height = 550))
      )
    )
  )
)

server <- function(input, output, session) {

  filtered <- reactive({
    d <- demo
    if (!identical(input$region, "All")) {
      d <- d[as.character(d$region) == input$region, , drop = FALSE]
    }
    d
  })

  output$map <- renderLeaflet({
    d <- filtered()
    vals <- d[[input$metric]]

    if (input$metric %in% c("difference", "pct_difference")) {
      lim <- max(abs(vals), na.rm = TRUE)
      pal <- colorNumeric(
        palette = "RdBu",
        domain = c(-lim, lim),
        reverse = TRUE,
        na.color = "transparent"
      )
    } else {
      pal <- colorNumeric(
        palette = "viridis",
        domain = vals,
        na.color = "transparent"
      )
    }

    label <- paste0(
      "<strong>", as.character(d$region), "</strong><br/>",
      "Model A: ", round(d$model_a, 1), "<br/>",
      "Model B: ", round(d$model_b, 1), "<br/>",
      "Difference: ", round(d$difference, 1), "<br/>",
      "Percent difference: ", round(d$pct_difference, 1), "%<br/>",
      "Buildings: ", d$buildings
    )

    leaflet(d) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addCircleMarkers(
        lng = ~lon,
        lat = ~lat,
        radius = input$radius,
        stroke = FALSE,
        fillOpacity = 0.85,
        fillColor = pal(vals),
        label = lapply(label, HTML)
      ) |>
      addLegend(
        "bottomright",
        pal = pal,
        values = vals,
        title = input$metric
      )
  })

  output$scatter <- renderPlot({
    d <- filtered()
    plot(
      d$model_a,
      d$model_b,
      xlab = "Model A posterior mean population",
      ylab = "Model B posterior mean population",
      main = "Cell-level comparison"
    )
    abline(0, 1, lty = 2)

    ok <- is.finite(d$model_a) & is.finite(d$model_b)
    if (sum(ok) > 2) {
      abline(lm(d$model_b[ok] ~ d$model_a[ok]), lty = 3)
    }
  })

  output$regional <- renderTable({
    x <- split(demo, demo$region)
    out <- lapply(names(x), function(r) {
      d <- x[[r]]
      a <- sum(d$model_a)
      b <- sum(d$model_b)
      data.frame(
        region = r,
        model_a_total = round(a, 1),
        model_b_total = round(b, 1),
        difference = round(b - a, 1),
        percent_difference = round(100 * (b - a) / pmax(a, 1e-8), 2)
      )
    })
    do.call(rbind, out)
  })

  output$diagnostics <- renderPlot({
    d <- filtered()

    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(2, 2))

    hist(
      d$difference,
      breaks = "FD",
      main = "B - A differences",
      xlab = "Population difference"
    )

    hist(
      d$pct_difference,
      breaks = "FD",
      main = "Percent differences",
      xlab = "Percent"
    )

    plot(
      d$model_a,
      d$difference,
      xlab = "Model A population",
      ylab = "B - A",
      main = "Difference vs magnitude"
    )
    abline(h = 0, lty = 2)

    plot(
      d$cv_a,
      d$cv_b,
      xlab = "Model A CV",
      ylab = "Model B CV",
      main = "Uncertainty comparison"
    )
    abline(0, 1, lty = 2)
  })
}

shinyApp(ui, server)
