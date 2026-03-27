Get_sequence <- function(AC_number) {
# Gets nucleotide sequences from Genbank Accession number
  if (is.na(AC_number) == TRUE) { NA } else {
    entrez_fetch(
      db      = "nuccore",
      id      = AC_number,
      rettype = "fasta",
      retmode = "text"
    ) |> 
      str_split_1("\n") |>
      _[-1] |>
      paste(collapse = "")
  }
}

write_to_fasta <- function(tib, path) {
# Writes two column tib as a .fa file
# Tibble structure must be :
# fasta name | sequence
  tib$accession <- paste0(">", tib$accession)
  c(rbind(tib$accession, tib$sequences)) |> write(file = path)
}

arrayspecs <- function (A, p, k, sep = NULL) {
  if (!is.matrix(A) && !is.data.frame(A)) 
    stop("A must be a data frame or matrix")
  dnames <- dimnames(A)
  n <- length(unlist(A))/(p * k)
  if (k < 2) 
    stop("One-dimensional data cannot be used")
  if (all(is.na(match(c(n, n * p), NROW(A))))) 
    stop("Matrix dimensions do not match input")
  specimens <- aperm(array(t(A), c(k, p, n)), c(2, 1, 3))
  dimnames(specimens)[[3]] <- dnames[[1]]
  col.names <- dnames[[2]]
  if (is.null(sep)) 
    col.names <- NULL
  if (!is.null(col.names)) {
    if (sep == "") {
      options(warn = -1)
      split.lab <- sort(unique(unlist(strsplit(col.names, 
                                               split = ""))))
      split.lab <- split.lab[length(split.lab)]
      a <- strsplit(col.names, split = split.lab)
      a <- a[sapply(a, length) == 2]
      a <- na.omit(as.numeric(unlist(a)))
      if (length(unlist(a)) == p) {
        split.no <- as.character(a[length(a)])
        b <- strsplit(col.names, split = split.no)
        no.char <- sapply(b, nchar)
        dim.tag <- which(no.char == min(no.char))
        b <- unlist(b[dim.tag])
        rn <- a
        cn <- b
      }
      else rn <- cn <- NULL
      options(warn = 0)
    }
    else {
      sep = paste("[", noquote(sep), "]")
      a <- strsplit(col.names, split = sep)
      if (length(unlist(a)) == 2 * k * p) {
        b <- simplify2array(a)
        rn <- unique(b[1, ])
        cn <- unique(b[2, ])
        if (length(rn) < length(cn)) {
          tmp <- cn
          cn <- rn
          rn <- tmp
        }
      }
      else rn <- cn <- NULL
    }
    dnames2 <- list(rn = rn, cn = cn)
  }
  if (!is.null(sep)) 
    dimnames(specimens)[1:2] <- dnames2
  return(specimens)
}

percentages <- function(pca) {
  pca |>
    (\(x) (x$sdev^2)/sum(x$sdev^2)*100)() |> # Compute percentage of variance drom sdev
    (\(x) tibble(Comps = 1:length(x), Values = x))() # Convert to tibble
}

base_breaks_x <- function(x) {
  b <- pretty(x)
  d <- data.frame(y = -Inf, yend = -Inf, x = min(b), xend = max(b))
  list(
    geom_segment(
      data        = d,
      mapping     = aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE
    ),
    scale_x_continuous(breaks = b)
  )
}

base_breaks_y <- function(x) {
  b <- pretty(x)
  d <- data.frame(x = -Inf, xend = -Inf, y = min(b), yend = max(b))
  list(
    geom_segment(
      data        = d,
      mapping     = aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE
    ),
    scale_y_continuous(breaks = b))
}

add_tree <- function(plot, tip_scores, tree, model = "BM") {
  
  data <- matrix(
    unlist(tip_scores |> _[, 2:3]),
    ncol = 2,
    dimnames = list(tip_scores$IDs)
  )
  
  ancestral_scores <- mvgls(pca_axes ~ 1, data = list("pca_axes" = data),
                            tree = tree, model = "BM", REML = TRUE, error = TRUE) |> 
    ancestral()
  colnames(ancestral_scores) <- c("x", "y")
  
  coordinates <- rbind(data[tree$tip.label, ], ancestral_scores) |>
    data.frame() |> 
    dplyr::rename(pcx = x, pcy = y) |> 
    dplyr::mutate(row = row_number())
  
  tree_edges <- tree$edge
  colnames(tree_edges) <- c("start", "end")
  
  pca_edges <- tree_edges |>
    as_tibble() |>
    left_join(coordinates, by = c("start" = "row")) |>
    dplyr::rename(x = pcx, y = pcy) |> 
    left_join(coordinates, by = c("end" = "row")) |>
    dplyr::rename(xend = pcx, yend = pcy)
  
  plot +
    geom_point(data = ancestral_scores |> data.frame(), aes(x = x, y = y), shape = 1) +
    geom_segment(data = pca_edges, aes(x = x, xend = xend, y = y, yend = yend))
}

PCA_tree_plot <- function(pca,
                          axis1 = 1, axis2 = 2,
                          colour_group = NULL, colour_name = NULL, palette = NULL,
                          labels       = TRUE,
                          phylo = NULL, model = "BM",
                          legend.position = "none") {
  
  tb <- tibble(
    IDs = rownames(pca$x),
    x   = pca$x[, axis1],
    y   = pca$x[, axis2]
  ) |>
    mutate(colour_group = colour_group)
  
  percent <- percentages(pca) |> pull(Values) |> round(2)
  
  p <- ggplot()
  
  if (!is.null(phylo)) {
    p <- add_tree(p, tip_scores = tb, tree = phylo, model = model)
  }
  
  p <- p +
    geom_point(
      data = tb,
      aes(x = x, y = y, colour = colour_group),
      size = 3
    )
  
  if (is.character(labels) == TRUE) {
    p <- p + geom_text(
      data = tb,
      aes(x = x, y = y, label = labels, colour = colour_group)
    )
  } else if (isTRUE(labels)) {
    p <- p + geom_text(
      data = tb,
      aes(x = x, y = y, label = IDs, colour = colour_group)
    )
  }
  
  p +
    scale_colour_discrete(name = colour_name, palette = palette) +
    base_breaks_x(tb$x) +
    base_breaks_y(tb$y) +
    geom_vline(
      xintercept = c(0),
      color      = "black",
      linetype   = "dashed",
      linewidth  = .12
    ) +
    geom_hline(
      yintercept = c(0),
      color      = "black",
      linetype   = "dashed",
      linewidth  = .12
    ) +
    coord_fixed(ratio = 1) +
    labs(
      x = paste0("PC", as.character(axis1), " = ", as.character(percent[[axis1]]), "%"),
      y = paste0("PC", as.character(axis2), " = ", as.character(percent[[axis2]]), "%")
    ) +
    theme(
      legend.position  = legend.position,
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}