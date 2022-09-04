\newcommand{\pkgfun}{\ifelse{\Sexpr[results=rd,stage=render]{length(find.package("#1", quiet = TRUE, verbose = FALSE)) > 0}}{\code{\link[#1:#2]{#1::#2()}}}{\code{#1::#2()}}}
