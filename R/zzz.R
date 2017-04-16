.onAttach <- function(libname, pkgname) {
    msg <- c(
        "\n*** Deprecation warning ***\n",
        "The ", pkgname, " package is deprecated and will be removed from\n",
        "Bioconductor 3.6.\n"
    )
    msg <- paste(msg, collapse="")
    .Deprecated(msg=msg)

    invisible(NULL)
}
