# Limit memory to 30GB:
if ("unix" %in% utils::installed.packages()[, 1]) {
    cat("Memory limit set to 30Gb\n\n")
    lims <- unix::rlimit_as(30*1024*1024*1024)
    rm(lims)

} else if ("ulimit" %in% utils::installed.packages()[, 1]) {
    cat("Memory limit set to 30Gb\n\n")
    lims <- ulimit::memory_limit(30*1024)
    rm(lims)

} else {
    warning("Memory usage will not be restricted. Install either 'unix' or 'ulimit' to limit memory use.")
}
