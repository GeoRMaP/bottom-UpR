library(shiny)
library(leaflet)

set.seed(2026)

make_demo_grid <- function() {
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

  ppb <- exp(
    1.0 +
      0.18 * sin(g$lon * 1.1) -
      0.14 * cos(g$lat * 0.9) +
      rnorm(n, 0, 0.10)
  )

  g$mean <- g$buildings * ppb
  g$sd <- pmax(0.5, 0.18 * g$mean + 0.3)
  g$lower <- pmax(0, g$mean - 1.96 * g$sd)
  g$upper <- g$mean + 1.96 * g$sd
  g$cv <- g$sd / pmax(g$mean, .Machine$double.eps)
  g$median <- pmax(0, g$mean - 0.08 * g$sd)

  # Structural zeros.
  zero <- g$buildings == 0
  g[zero, c("mean", "sd", "lower", "median", "upper", "cv")] <- 0

  g
}

demo <- make_demo_grid()

ui <- fluidPage(
  titlePanel("bottom-UpR: Gridded Posterior Results Explorer"),
  sidebarLayout(
    sidebarPanel(
      helpText(
        "This bundled demo uses synthetic gridded posterior results. ",
        "It illustrates how fine-grid population outputs can be explored interactively."
      ),
      selectInput(
        "metric",
        "Map variable",
        choices = c(
          "Posterior mean population" = "mean",
          "Posterior median population" = "median",
          "Lower 95% interval" = "lower",
          "Upper 95% interval" = "upper",
          "Posterior SD" = "sd",
          "Coefficient of variation" = "cv",
          "Mapped buildings" = "buildings"
        ),
        selected = "mean"
      ),
      selectInput(
        "region",
        "Region",
        choices = c("All", levels(demo$region)),
        selected = "All"
      ),
      sliderInput(
        "radius",
        "Grid-cell display radius",
        min = 2,
        max = 8,
        value = 4,
        step = 1
      ),
      sliderInput(
        "opacity",
        "Cell opacity",
        min = 0.2,
        max = 1,
        value = 0.8,
        step = 0.1
      ),
      checkboxInput(
        "structural",
        "Highlight structural-zero cells",
        value = TRUE
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Map",
          leafletOutput("map", height = 650)
        ),
        tabPanel(
          "Distribution",
          plotOutput("hist", height = 500)
        ),
        tabPanel(
          "Summary",
          verbatimTextOutput("summary"),
          tableOutput("region_table")
        )
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
    metric <- input$metric
    vals <- d[[metric]]

    pal <- colorNumeric(
      palette = "viridis",
      domain = vals,
      na.color = "transparent"
    )

    labels <- paste0(
      "<strong>", as.character(d$region), "</strong><br/>",
      "Mean population: ", round(d$mean, 1), "<br/>",
      "95% interval: ", round(d$lower, 1), " - ", round(d$upper, 1), "<br/>",
      "CV: ", round(d$cv, 3), "<br/>",
      "Buildings: ", d$buildings
    )

    m <- leaflet(d) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addCircleMarkers(
        lng = ~lon,
        lat = ~lat,
        radius = input$radius,
        stroke = FALSE,
        fillOpacity = input$opacity,
        fillColor = pal(vals),
        label = lapply(labels, HTML)
      ) |>
      addLegend(
        "bottomright",
        pal = pal,
        values = vals,
        title = names(which(
          c(
            mean = "Posterior mean population",
            median = "Posterior median population",
            lower = "Lower 95% interval",
            upper = "Upper 95% interval",
            sd = "Posterior SD",
            cv = "Coefficient of variation",
            buildings = "Mapped buildings"
          ) == switch(
            metric,
            mean = "Posterior mean population",
            median = "Posterior median population",
            lower = "Lower 95% interval",
            upper = "Upper 95% interval",
            sd = "Posterior SD",
            cv = "Coefficient of variation",
            buildings = "Mapped buildings"
          )
        ))
      )

    if (isTRUE(input$structural)) {
      z <- d[d$buildings == 0, , drop = FALSE]
      if (nrow(z)) {
        m <- m |>
          addCircleMarkers(
            data = z,
            lng = ~lon,
            lat = ~lat,
            radius = input$radius + 1,
            stroke = TRUE,
            weight = 1.5,
            color = "black",
            fill = FALSE,
            label = "Structural zero: no mapped buildings"
          )
      }
    }

    m
  })

  output$hist <- renderPlot({
    d <- filtered()
    x <- d[[input$metric]]

    hist(
      x,
      breaks = "FD",
      main = paste("Distribution of", input$metric),
      xlab = input$metric
    )
    abline(v = mean(x, na.rm = TRUE), lty = 2)
  })

  output$summary <- renderPrint({
    d <- filtered()

    cat("Grid cells:", nrow(d), "\n")
    cat("Structural-zero cells:", sum(d$buildings == 0), "\n")
    cat("Posterior mean total population:", round(sum(d$mean), 1), "\n")
    cat("Mean cell CV:", round(mean(d$cv[d$buildings > 0], na.rm = TRUE), 3), "\n\n")

    print(summary(d[[input$metric]]))
  })

  output$region_table <- renderTable({
    x <- split(demo, demo$region)
    out <- lapply(names(x), function(r) {
      d <- x[[r]]
      data.frame(
        region = r,
        cells = nrow(d),
        buildings = sum(d$buildings),
        posterior_mean_population = round(sum(d$mean), 1),
        mean_cv = round(mean(d$cv[d$buildings > 0], na.rm = TRUE), 3)
      )
    })
    do.call(rbind, out)
  })
}

shinyApp(ui, server)
