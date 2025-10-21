
# скрипт для установки необходимых пакетов

packages <- readLines("requirements.txt")

install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(packages, install_if_missing))
