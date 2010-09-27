#.onLoad <- function(libname, pkgname) {
#  library.dynam(pkgname)
#}

.onUnload <- function(libpath)
{
    library.dynam.unload("ppstat", libpath)
}
