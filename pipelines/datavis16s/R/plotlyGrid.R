#' @title Add Plotly data export to Plotly graph
#'
#' @description All functions create an output html plot  with link which sends the data to a grid in the plotly chart
#'  studio.
#'
#'  \code{plotlyGrid} takes in a ggplot or plotly object and creates an output html plotly plot.
#'
#' @param pplot plotly or ggplot object
#' @param filename output filename (fullpath)
#' @param data data frame to export to plotly grid (optional for plotlyGrid)
#' @param title title of html page
#' @param outlib (Optional) name of external lib directory for non-selfcontained html.
#' Useful for multiple graphs sharing the same lib.
#'
#' @importFrom htmlwidgets prependContent saveWidget
#' @importFrom plotly plotly_data plotly_build ggplotly
#'
#' @return html plot is saved to filename. external libraries are saved to outlib in same directory as filename.
#' Invisibly returns the plotly html widget.
#'
#' @source [plotlyGrid.R](../R/plotlyGrid.R)
#' @rdname plotlyGrid
#' @name plotlyGrid
#'
plotlyGrid <- function(pplot, filename, data=NULL, title=NULL, outlib="lib") {

  if ("ggplot" %in% class(pplot)) {

    withCallingHandlers({
      pplot <- ggplotly(pplot)
    }, message=function(c) {
      if (startsWith(conditionMessage(c), "We recommend that you use the dev version of ggplot2"))
        invokeRestart("muffleMessage")
    })
  }
  pp <- plotly_build(pplot)

  if (is.null(data)) {
    data <- plotly_data(pp)
  }

  list2env(gridCode(data), envir=environment())

  pp <- htmlwidgets::appendContent(pp, html)
  pp <- htmlwidgets::appendContent(pp,javascript)

  if (is.null(title)) {
    title <- tools::file_path_sans_ext(basename(filename))

  }

  outfile <- file.path(tools::file_path_as_absolute(dirname(filename)), basename(filename))
  outlib <- file.path(dirname(outfile), basename(outlib))
  logoutput(paste("Saving plot to", outfile))
  plotlywidget <- plotly::config(pp,  cloud=T, edits = list(titleText=T, legendText=T, legendPosition=T, axisTitleText=T))
  saveWidget(plotlywidget, file=outfile , selfcontained = FALSE, title=title, libdir=outlib)
  invisible(plotlywidget)
}

#' @title Add Plotly data export to plain html
#'
#' @description \code{htmlGrid} takes in an html tag object.
#'
#' @param ht html tagList
#' @param jquery should we load jquery
#' @param styletags html object with style tags for the tagList.
#'
#' @importFrom htmltools tagList tags
#' @importFrom rmarkdown html_dependency_jquery
#'
#' @details If jquery is needed, we use jquery-1.11.3 from the rmarkdown library.  We also use
#'  shiny's bootstrap-3.3.7 css to style the text elements.
#'
#' @rdname plotlyGrid
#'
htmlGrid <- function(ht, filename, data, jquery = FALSE, title=NULL, outlib="lib", styletags=NULL) {

  list2env(gridCode(data), envir=environment())
  tl <- tagList(javascript, ht)

  if (jquery) {
    jq <- html_dependency_jquery()
    tl <- htmltools::attachDependencies(tl, jq, append=TRUE)
  }

  mc <- shiny::bootstrapLib()
  mc <- htmltools::htmlDependency("bootstrap-css", version=mc$version, src=mc$src, stylesheet = mc$stylesheet, all_files = F)
#  mc <- htmltools::htmlDependency("bootstrap-css", "3.3.7", c(file=file.path(find.package("shiny"), "www/shared/bootstrap/css")), stylesheet = "bootstrap.min.css", all_files = F)
  tl <- htmltools::attachDependencies(tl, mc, append=TRUE)
  if (!is.null(title)) {
    tl <- tags$div(class="container-fluid", style="max-width:1200px", tags$h2(title), tags$p(html), tl)
  } else {
    tl <- tags$div(class="container-fluid", tags$p(html), tl)
  }

  if (!is.null(styletags)) {
    tl <- tagList(styletags, tl)
  }

  outfile <- file.path(tools::file_path_as_absolute(dirname(filename)), basename(filename))
  outlib <- file.path(dirname(outfile), basename(outlib))

  logoutput(paste("Saving plot to", outfile))
  htmltools::save_html(tl, file= outfile, lib=outlib)
}


#' @title Format plotly grid code
#'
#' @description Format data according to here: \url{https://plot.ly/export/}
#'
#' @param data data to populate plotly grid
#'
#' @return list of 2 values:
#' \describe{
#'   \item{html}{html for plotly export link}
#'   \item{javascript}{js function for exporting data}
#' }
#'
#' @importFrom htmltools HTML
#'
#' @source [plotlyGrid.R](../R/plotlyGrid.R)
#'
gridCode <- function(data) {

  ll <- as.list(data)

  nn <- names(ll)
  ll <- lapply(nn, function(x) { mm <- ll[[x]]; ll[[x]] <- NULL; return(list(data=mm));  } )
  names(ll) <- nn
  ll <- jsonlite::toJSON(ll)

  html <- HTML(text = '<a href=\"#\" id=\"plotly-data-export\" target=\"_blank\" style=\"font-family:\'Open Sans\',sans-serif;\">Export Data to Plotly</a>')

  javascript <- HTML(paste('<script>',
                           'function getPlotlyGridData(){',
                           '  return {',
                           paste0('cols: ', ll),
                           '  }',
                           '}',
                           "$(\'#plotly-data-export\').click(function(){",
                           "  var hiddenForm = $(\'<div id=\"hiddenform\" \'+",
                           "                       \'style=\"display:none;\">\'+",
                           "                       \'<form action=\"https://chart-studio.plotly.com/datagrid\" \'+",
                           "                       \'method=\"post\" target=\"_blank\">\'+",
                           "                       \'<input type=\"text\" \'+",
                           "                       \'name=\"data\" /></form></div>\')",
                           "  .appendTo(\'body\');",
                           "  var dataGrid = JSON.stringify(getPlotlyGridData())",
                           "  hiddenForm.find(\'input\').val(dataGrid);",
                           "  hiddenForm.find(\'form\').submit();",
                           "  hiddenForm.remove();",
                           "  return false;",
                           "});",
                           "</script>", sep = "\n"))
  return(list(html=html, javascript=javascript))
}


#' Save an HTML object to a file
#'
#' @param html HTML content to print
#' @param file File to write content to
#' @param background Background color for web page
#' @param libdir Directory to copy dependencies to
#' @param bodystyle html style string
#'
#' @return save html to file
#'
#' @source [plotlyGrid.R](../R/plotlyGrid.R)
#'
save_fillhtml <- function (html, file, background = "white", libdir = "lib", bodystyle="")
{
  dir <- dirname(file)
  oldwd <- setwd(dir)
  on.exit(setwd(oldwd), add = TRUE)
  rendered <- htmltools::renderTags(html)
  deps <- lapply(rendered$dependencies, function(dep) {
    dep <- htmltools::copyDependencyToDir(dep, libdir, FALSE)
    dep <- htmltools::makeDependencyRelative(dep, dir, FALSE)
    dep
  })
  html <- c("<!DOCTYPE html>", "<html>", "<head>", "<meta charset=\"utf-8\"/>",
            htmltools::renderDependencies(deps, c("href", "file")), rendered$head,
            "</head>", sprintf("<body style=\"background-color:%s;%s\">",
                               htmltools::htmlEscape(background), bodystyle), rendered$html, "</body>",
            "</html>")
  writeLines(html, file, useBytes = TRUE)
}

