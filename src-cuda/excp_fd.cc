double excp_fd (double pfhs, double dens) {

  double fd = (4.0 - 2.0 * pfhs) / ((1.0 - pfhs) * (1.0 - pfhs) * (1 - pfhs))
    * (pfhs / dens);
  return fd;
}
