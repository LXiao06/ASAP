#include <Rcpp.h>

#ifdef _WIN32
// Windows: no-op (fontconfig warnings don't occur)
#else
#include <unistd.h>
#include <fcntl.h>
#endif

//' Suppress C-level stderr output
//'
//' Redirects file descriptor 2 (C stderr) to /dev/null and returns
//' the saved file descriptor so it can be restored later.
//' This suppresses C-library warnings (e.g. fontconfig's
//' "using without calling FcInit()") that R's sink() cannot intercept.
//'
//' @return Integer file descriptor to pass to \code{restore_stderr()},
//'   or -1 on Windows where this is a no-op.
//' @keywords internal
// [[Rcpp::export]]
int suppress_stderr() {
#ifndef _WIN32
  int saved_fd = dup(STDERR_FILENO);
  int devnull = open("/dev/null", O_WRONLY);
  if (devnull >= 0) {
    dup2(devnull, STDERR_FILENO);
    close(devnull);
  }
  return saved_fd;
#else
  return -1;
#endif
}

//' Restore C-level stderr output
//'
//' Restores file descriptor 2 from a previously saved descriptor
//' returned by \code{suppress_stderr()}.
//'
//' @param saved_fd Integer file descriptor from \code{suppress_stderr()}.
//' @keywords internal
// [[Rcpp::export]]
void restore_stderr(int saved_fd) {
#ifndef _WIN32
  if (saved_fd >= 0) {
    dup2(saved_fd, STDERR_FILENO);
    close(saved_fd);
  }
#endif
}
