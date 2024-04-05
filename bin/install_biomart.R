options(install.packages.compile.from.source = "always")
.libPaths("/workspace/Goatee/lib")
install.packages("BiocManager", lib="/workspace/Goatee/lib")
BiocManager::install("biomaRt")