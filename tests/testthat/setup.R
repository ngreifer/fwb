#Open a null graphics device for the whole test run.
#
#`plot.fwb()` draws to the active device. With none open, R starts the default one,
#which writes an `Rplots.pdf` into tests/testthat/. Opening the null device here means
#nothing is ever written.
#
#No teardown is registered: the device closes when the test process exits, and tests
#that draw use `local_null_device()` (see helpers.R), which opens and closes one
#around themselves.
grDevices::pdf(NULL)
