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
#' @importFrom htmlwidgets prependContent appendContent saveWidget
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
  pp <- htmlwidgets::prependContent(pp, html)
  pp <- htmlwidgets::appendContent(pp,javascript)

  if (is.null(title)) {
    title <- tools::file_path_sans_ext(basename(filename))

  }

  outfile <- file.path(tools::file_path_as_absolute(dirname(filename)), basename(filename))
  outlib <- file.path(dirname(outfile), basename(outlib))
  logoutput(paste("Saving plot to", outfile))
  pp <- plotly::config(pp,  cloud=T, edits = list(titleText=T, legendText=T, legendPosition=T, axisTitleText=T))
  pp$sizingPolicy$padding <- 40
  saveWidget(pp, file=outfile , selfcontained = FALSE, title=title, libdir=outlib)
  invisible(pp)
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

  html <- HTML(text = '<div id=\"plotly-data-link\" style=\"width: 100%; height: 40px;\"> <a href=\"#\" id=\"plotly-data-export\" target=\"_blank\" style=\"font-family:\'Open Sans\',sans-serif;\">Export Data to Plotly</a> </div>')


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
    dep$all_files <- FALSE
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

