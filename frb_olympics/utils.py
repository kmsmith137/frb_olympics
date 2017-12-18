# Hmm, I'll probably put this somewhere else eventually...

# Note: cut-and-paste from "simpulse"
def dispersion_delay(dm, freq_MHz):
    """Returns dispersion delay in seconds, for given DM and frequency."""
    return 4.148806e3 * dm / (freq_MHz * freq_MHz);
