\newcommand{\pkgfun}{\Sexpr[results=rd,stage=render]{if (rlang::is_installed("#1")) r"(\code{\link[#1:#2]{#1::#2()}})" else r"(\code{#1::#2()})"}}
